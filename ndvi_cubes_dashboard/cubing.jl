using GLMakie, YAXArrays
using Zarr, DimensionalData
using Colors
using TileProviders
using Tyler
using Extents
using MapTiles
function to_web_mercator(lo,lat)
    return Point2f(MapTiles.project((lo,lat), MapTiles.wgs84, MapTiles.web_mercator))
end

function getStats(target, prediction)
    outstats = zeros(size(target)[2:end]..., 9) # 9 is the amount of metrics from nseparts
    for i in axes(target, 2), j in axes(target, 3)
        yₘ, ŷₘ = target[:, i, j], prediction[:, i, j]
        outstats[i, j, :] .= nseparts(yₘ, ŷₘ)
    end
    #Δ□tp = target .- prediction 
    return outstats
end

"""
nseparts(yₘ, ŷₘ)
Compute NSE and components
"""
function nseparts(yₘ, ŷₘ)
    μ₀, μ̂ = mean(skipmissing(yₘ)), mean(skipmissing(ŷₘ))
    σ₀, σ̂ = sqrt(var(skipmissing(yₘ))), sqrt(var(skipmissing(ŷₘ)))
    idxs = (.!ismissing.(yₘ)) .& (.!ismissing.(ŷₘ))
    r = sum(idxs) > 0 ? cor(yₘ[idxs], ŷₘ[idxs]) : NaN
    α = σ̂ / σ₀
    β = (μ̂ - μ₀) / σ₀
    cov_val = sum(idxs) > 0 ? cov(yₘ[idxs], ŷₘ[idxs]) : NaN
    return [2 * α * r - α^2 - β^2, α, r, β, μ₀, σ₀, μ̂, σ̂, cov_val]
end



ds = open_dataset(joinpath(@__DIR__, "cubed_ndvi_stacked.zarr"))
c = ds["data"]
lon_lat = ds["lon_lat"]
lons = lon_lat.data[1,64,:]
lats = lon_lat.data[2,64,:]
locs_cubes = Point2f.(lons, lats)
cubes_names = lookup(c, :cube_name)

wm_points = to_web_mercator.(lons,lats)
init_index = findall(x->x=="28PCU0804", cubes_names)[1]
# generate point instances
ndvi_points = [Point3f(i,j,k) for i in range(1,32,10) for j in 1:128 for k in 1:128];

# Observables
name_var = Observable("ndvi_target")
idx_pos = Observable(init_index)

mv_point = @lift(wm_points[$idx_pos])
# lift for cube name !
c_name = @lift(cubes_names[$idx_pos])  #Observable("28PCU0804")
#ndvi_data = @lift(c[ndvi=At($name_var)]);
#mn_cube = @lift($ndvi_data[cube_name=At($c_name)].data[:,:,end:-1:1]);

mn_cube = lift(name_var, c_name) do var_name, cube_n
    ndvi_data = c[ndvi=At(var_name)]
    ndvi_data[cube_name=At(cube_n)].data[:,:,end:-1:1]
end

ndvi_values = @lift([$mn_cube[i,j,k] for i in 1:10 for j in 1:128 for k in 1:128]);

cmap = :Spectral
fs=1
ilon = minimum(lons) - 1
flon = maximum(lons) + 2
ilat = minimum(lats) - 1
flat = maximum(lats) + 0.2

extent = Extent(X = (ilon, flon), Y = (ilat, flat));

with_theme(theme_minimal()) do
    fig = Figure(; figure_padding=(5,5,5,0), size = (9*60,16*60), backgroundcolor="#f5f5f3")
    menu = Menu(fig, options = ["target", "prediction"])
    #menu_names = Menu(fig, options = cubes_names)
    g_lay = GridLayout(fig[2,1])
    ax_inset = Axis(g_lay[1:2,1])

    ax_series = Axis(g_lay[2,2:5]; yaxisposition=:right)
    #hideydecorations!(ax_series)
    hidespines!.(ax_series, :l)
    axs = [Axis(g_lay[1,i], aspect=1) for i in 2:5]
    hidedecorations!.(axs)
    hidespines!.(axs)
    [heatmap!(axs[k], rand(128,128)) for k in 1:4]

    colsize!(g_lay, 1, Auto(2))
    rowsize!(g_lay, 1, Auto(1.5))

    colgap!(g_lay,5)
    rowgap!(g_lay,5)

    m = Tyler.Map(extent; provider=CartoDB(), figure=fig,
        axis=ax_inset, scale=3)

    p = scatter!(m.axis, wm_points; marker=:rect, markersize=10, color=(:tan1, 0.85),
        strokewidth=0.75, strokecolor=:black)

    scatter!(m.axis, mv_point; marker=:rect, markersize=16, color=(:white, 0.1),
        strokewidth=2, strokecolor=:dodgerblue)

    hidedecorations!(ax_inset)
    hidespines!(ax_inset)
    ax = LScene(fig[1,1]; show_axis=false)
    plt = meshscatter!(ndvi_points; marker=Rect3f(Vec3f(-0.5),Vec3f(1)),
        color=ndvi_values, colormap=cmap, shading=FastShading, colorrange=(0,1),
        lowclip=:grey9, nan_color= (:grey10, 0.05),
        markersize=Vec3f(0.8*3,0.99,0.99))

    Colorbar(fig[1,1, Bottom()], plt; vertical=false,
        label=@lift(rich("NDVI  ", rich("$($c_name)", font=:bold, color=:tan1))),
        width=Relative(0.3), tellwidth=false, halign=0)

    # fig[1, 1] = vgrid!(Label(fig, "", width = nothing), menu_names;
    #     tellheight = false,tellwidth = false, width=150, halign=1, valign=1)

    fig[1, 1] = vgrid!(Label(fig, "", width = nothing), menu;
        tellheight = false,tellwidth = false, width=150, halign=1, valign=1)
    colgap!(fig.layout, 5)
    rowgap!(fig.layout, 2)
    rowsize!(fig.layout, 2, Auto(0.3))

    on(menu.selection) do s
        new_s = s == "target" ? "ndvi_target" : "ndvi_pred"
        name_var[] = new_s
    end
    # events && interactions
    # on(menu_names.selection) do s_name
    #     c_name[] = s_name
    # end
    on(events(m.axis).mousebutton, priority = 2) do event
        if event.button == Mouse.left && event.action == Mouse.press
            # move marker
            plt, i = Makie.pick(m.axis)
            if plt == p
                idx_pos[] = i
                notify(idx_pos)
                notify(name_var)
                return Consume(true)
            end
            #center!(ax.scene)
        end
        return Consume(false)
    end
    fig
end

# using YAXArrays, Zarr, DimensionalData
# using GLMakie
# GLMakie.activate!(float=true)

# ds = open_dataset("/Users/lalonso/Documents/LittlePieces/tmp_ndvi_cubes/cubed_ndvi_stacked.zarr")
# c = ds["data"]
# cubes_names = lookup(c, :cube_name)
# # generate point instances
# ndvi_points = [Point3f(i,j,k) for i in range(1,32,10) for j in 1:128 for k in 1:128];

# sample_data=[]
# ids = ["28PCU0804", "28PCU0305", "28PCV1209", "30STC0601", "36RTU2305", "34HCH0705"]
# for c_name in ids
#     ndvi_t = c[ndvi=At("ndvi_target")]
#     ndvi_p = c[ndvi=At("ndvi_pred")]
#     new_d = ndvi_t[cube_name=At(c_name)].data[:,:,end:-1:1]
#     new_p = ndvi_p[cube_name=At(c_name)].data[:,:,end:-1:1]
#     push!(sample_data, [new_d, new_p])
# end
# sample_data =[s for s in sample_data];
# jldsave(joinpath(@__DIR__, "sample_cubes.jld2"); target_prediction = sample_data)

# data_jld = jldopen(joinpath(@__DIR__, "sample_cubes.jld2"), "r")
# data_array = data_jld["target_prediction"]
