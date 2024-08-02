using GLMakie, YAXArrays
using Zarr, DimensionalData
using Colors
using TileProviders
using Tyler
using Extents
using MapTiles
using Statistics
using NaNStatistics
using Dates
function to_web_mercator(lo,lat)
    return Point2f(MapTiles.project((lo,lat), MapTiles.wgs84, MapTiles.web_mercator))
end

function getStats(target, prediction)
    outstats = zeros(size(target)[1:2]..., 10) # 9 is the amount of metrics from nseparts + rmse
    for i in axes(target, 1), j in axes(target, 2)
        yₘ, ŷₘ = target[i, j, :], prediction[i, j, :]
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
    rmse = sqrt(nanmean((yₘ .- ŷₘ).^2))
    yₘ = replace(yₘ, NaN=>missing)
    ŷₘ = replace(ŷₘ, NaN=>missing)

    μ₀, μ̂ = mean(skipmissing(yₘ)), mean(skipmissing(ŷₘ))
    σ₀, σ̂ = sqrt(var(skipmissing(yₘ))), sqrt(var(skipmissing(ŷₘ)))
    idxs = (.!ismissing.(yₘ)) .& (.!ismissing.(ŷₘ))
    r = sum(idxs) > 0 ? cor(yₘ[idxs], ŷₘ[idxs]) : NaN
    α = σ̂ / σ₀
    β = (μ̂ - μ₀) / σ₀
    cov_val = sum(idxs) > 0 ? cov(yₘ[idxs], ŷₘ[idxs]) : NaN
    return [2 * α * r - α^2 - β^2, α, r, β, μ₀, σ₀, μ̂, σ̂, cov_val, rmse]
end
function linkscenes!(connected_axes)
    scenes = [ax.scene for ax in connected_axes]
    locks = [Observable(false) for _ in scenes]
    for i in eachindex(scenes)
        cam = scenes[i].camera_controls
        onany(cam.eyeposition, cam.lookat, cam.upvector) do eyepos, lookat, up
            locks[i][] && return
            for j in eachindex(scenes)
                i == j && continue
                locks[j][] = true
                update_cam!(scenes[j], eyepos, lookat, up)
                locks[j][] = false
            end
        end
    end
end
# add rmse and remove alpha
# add msc per slice in time series plot

ds = open_dataset(joinpath(@__DIR__, "ndvi_clean_extremes_4534_cubes.zarr/"))

tempo = ds["tempo_target"]
extremes = ds["extremes_10"]

c = ds["ndvi"]


lon_lat = ds["lon_lat"]
lons = lon_lat.data[1,64,:]
lats = lon_lat.data[2,64,:]
locs_cubes = Point2f.(lons, lats)
cubes_names = lookup(c, :cube_name)

wm_points = to_web_mercator.(lons,lats)
init_index = findall(x->x=="36KVF0050", cubes_names)[1]

# generate point instances
ndvi_points = [Point3f(i,j,k) for i in range(1,32,10) for j in 1:128 for k in 1:128];

# Observables
time_scroll = Observable(0)

idx_pos = Observable(init_index)
mv_point = @lift(wm_points[$idx_pos])
# lift for cube name !
c_name = @lift(cubes_names[$idx_pos])


mn_cube = lift(c_name) do cube_n
    ndvi_data = c[variable=At("ndvi_target")]
    ndvi_data[cube_name=At(cube_n)].data[:,end:-1:1,:]
end

mn_cube_p = lift(c_name) do cube_n
    ndvi_data = c[variable=At("ndvi_pred")]
    ndvi_data[cube_name=At(cube_n)].data[:,end:-1:1,:]
end

d_str_o = lift(c_name) do cube_n
    t_xticks = tempo[cube_name=At(cube_n)].data[:]
    t_str = string.(t_xticks)
    string.(Date.(t_str, DateFormat("yyyymmdd")))
end

extremes_10  = lift(c_name) do cube_n
    ex_factor = extremes[cube_name=At(cube_n)].data[:]
    replace(x->x>0 ? 1 : x, ex_factor)
end

μ_target = @lift([nanmean($mn_cube[:,:,k]) for k in 1:10])
μ_pred = @lift([nanmean($mn_cube_p[:,:,k]) for k in 1:10])

diff_t_p = @lift(abs.($μ_target .- $μ_pred))

#msc
msc_cube = lift(c_name) do cube_n
    ndvi_data = c[variable=At("msc")]
    ndvi_data[cube_name=At(cube_n)].data[:,end:-1:1,:]
end
μ_msc = @lift([nanmean($msc_cube[:,:,k]) for k in 1:10])


ndvi_values = @lift([$mn_cube[j,k,i] for i in 1:10 for j in 1:128 for k in 1:128]);
ndvi_values_pred = @lift([$mn_cube_p[j,k,i] for i in 1:10 for j in 1:128 for k in 1:128]);

stats = @lift(getStats($mn_cube, $mn_cube_p))

nse = @lift($stats[:,:,1])
alpha = @lift($stats[:,:,2])
r_corr = @lift($stats[:,:,3])
beta = @lift($stats[:,:,4])
rmse_val = @lift($stats[:,:,10])

# scroll time slices
scrolled_ndvi_values = @lift(circshift($ndvi_values, -9*($time_scroll*128*128)))
scrolled_ndvi_values_pred = @lift(circshift($ndvi_values_pred, -9*($time_scroll*128*128)))

cmap = :Spectral_11 #resample_cmap(:gist_earth, 100)[20:end]
fs=1
ilon = minimum(lons) - 1
flon = maximum(lons) + 2
ilat = minimum(lats) - 1
flat = maximum(lats) + 0.2

extent = Extent(X = (ilon, flon), Y = (ilat, flat));

nan_colors = [(:grey10, 0.025), (:dodgerblue, 0.025), (:orangered, 0.025), (:transparent, 0.0)]
labels_colors = ["grey10", "blue", "reds", "transparent"]

nan_color_o = Observable(nan_colors[1])

#set_theme!(theme_dark())
function plot_cubes(; provider=CartoDB(:DarkMatter))
    fig = Figure(; figure_padding=(15,15,5,0), size = (1920÷1.5,1080÷1.5),
        #backgroundcolor=1.15colorant"gainsboro",
        fontsize=16,
        fonts = (; regular="CMU Serif")
        )
    low_lay = GridLayout(fig[2,1:3])

    g_lay2 = GridLayout(low_lay[1,1:3])
    g_lay = GridLayout(low_lay[1,4:5])

    ax_inset = Axis(g_lay2[1,1:3], ylabel = "NDVI", #yaxisposition=:right
        )

    ax_map = Axis(fig[1,1])
    scatterlines!(ax_inset, μ_target, label = "spatial mean target",
        markersize=15,linestyle=:dash, color=:tan1
        )
    scatterlines!(ax_inset, μ_pred, label = "spatial mean prediction",
        marker=:rect, markersize=15, linestyle=:dash, color=0.7*colorant"steelblue1",
        )
    scatterlines!(ax_inset, μ_msc, label = "msc",
        marker=:utriangle, markersize=15, linestyle=:dash,
        markercolor=:transparent, strokewidth=1,
        strokecolor=colorant"purple2",
        color = colorant"purple2",
        )
    barplot!(ax_inset, 1:10, diff_t_p, width=0.2, color=:brown1, 
        label = "abs(target - prediction)")

    Legend(g_lay2[1,1:3], ax_inset, nbanks=4, tellwidth=false,tellheight = false,
        halign=0.5, valign=1.2)

    ax_inset.xticklabelrotation = π / 4
    ax_inset.xticklabelalign = (:right, :center)
    lift(d_str_o) do n_str
        ax_inset.xticks = (1:10, n_str)
    end
    vlines!(ax_inset, 1:10; color = extremes_10, colormap=[:transparent, :red],
        colorrange=(0, 1))

    ylims!(ax_inset, 0, 1)
    xlims!(ax_inset, 0.5, 10.5)

    axs = [Axis(g_lay[1,i], aspect=1) for i in 1:4]

    hidedecorations!.(axs)
    hidespines!.(axs)
    hm1 = heatmap!(axs[1], nse, colorrange=(-0.5,1), colormap=:CMRmap, lowclip=:black)
    hm2 = heatmap!(axs[2], r_corr, colorrange=(-1,1), colormap=:Hiroshige,)
    hm3 = heatmap!(axs[3], rmse_val, colorrange=(0,0.3), colormap=:thermal, highclip=:gold)
    hm4 = heatmap!(axs[4], beta, colorrange=(-2.5,2.5), colormap=:Spectral_11,
        lowclip=:red4, highclip=:midnightblue)
    linkaxes!.(axs[1:4]...)

    Colorbar(g_lay[1,1], hm1; vertical=false,
        label=L"NSE",
        width=Relative(0.7), tellwidth=false,tellheight = false,
        halign=0.5,
        valign=1.05
        )
    cb = Colorbar(g_lay[1,2], hm2; vertical=false,
        label="r",
        width=Relative(0.7), tellwidth=false, tellheight = false,
        halign=0.5,
        valign=1.05
        )
    cb.ticks = [-1,0,1]
    Colorbar(g_lay[1,3], hm3; vertical=false,
        label=L"RMSE",
        width=Relative(0.7), tellwidth=false,tellheight = false,
        halign=0.5,
        valign=1.05
        )
    Colorbar(g_lay[1,4], hm4; vertical=false,
        label=L"\beta",
        width=Relative(0.7), tellwidth=false,tellheight = false,
        halign=0.5,
        valign=1.05
        )

    colgap!(g_lay,5)
    rowgap!(g_lay,5)

    m = Tyler.Map(extent; provider=provider, figure=fig, # :DarkMatter
        axis=ax_map, scale=3)

    p = scatter!(m.axis, wm_points; marker=:rect,
        markerspace=:data,
        markersize= 2560,
        color=(:dodgerblue, 0.5),
        strokewidth=0.5, strokecolor=(:black, 0.5))

    scatter!(m.axis, mv_point; marker=:rect, markersize=16, color=(:grey45, 0.1),
        strokewidth=2.5, strokecolor=:red)

    hidedecorations!(ax_map)
    hidespines!(ax_map)

    ax = LScene(fig[1,2]; show_axis=false)

    plt = meshscatter!(ndvi_points; marker=Rect3f(Vec3f(-0.5),Vec3f(1)),
        color=scrolled_ndvi_values, colormap=cmap, shading=FastShading, colorrange=(0,1),
        lowclip=:grey9, nan_color= nan_color_o,
        markersize=Vec3f(0.35*3,0.99,0.99))

    ax_p = LScene(fig[1,3]; show_axis=false)

    plt_p = meshscatter!(ax_p, ndvi_points; marker=Rect3f(Vec3f(-0.5),Vec3f(1)),
        color=scrolled_ndvi_values_pred, colormap=cmap, shading=FastShading, colorrange=(0,1),
        lowclip=:grey9, nan_color= nan_color_o,
        markersize=Vec3f(0.35*3,0.99,0.99))
    linkscenes!([ax, ax_p])

    Colorbar(fig[1,2, Bottom()], plt; vertical=false,
        label=@lift(rich("NDVI  ", rich("$($c_name)", font=:bold, color=:tan1))),
        width=Relative(0.3), tellwidth=false, halign=0)

    sg = SliderGrid(
            fig[1, 2, Bottom()],
            (label = "Time Scroller", range = 0:10, startvalue = 0, color_active=:grey55,
            color_inactive=(:grey55, 0.2), color_active_dimmed=:gold1,),
            width=Relative(0.5),
            halign=0.75,
            tellheight = false)
    connect!(time_scroll, sg.sliders[1].value)

    Label(fig[0,1:end], rich(rich("MAX PLANCK INSTITUTE\n", color=:grey55, font=:bold),
        rich("FOR BIOGEOCHEMISTRY", color=:grey45, fontsize=12*fs));
        tellwidth=false, halign=0.05,  justification=:right, padding=(0,0,10,10)
        )

    nan_menu = Menu(fig, options = zip(labels_colors, nan_colors),
        selection_cell_color_inactive=(:grey55, 0.2),
        cell_color_active=colorant"gold1",
        cell_color_hover=colorant"grey85",
        textcolor=:grey55,
    )
    fig[0,3] = vgrid!(
        Label(fig, "nan_color", color=:orange, width = nothing),
        nan_menu;
        tellheight = false, tellwidth=false, width = 100,
        halign=1, valign=1.1)

    Label(g_lay[1,1:4, Bottom()], rich(rich("Visualization by ", color=:grey45), rich("Lazaro Alonso\n ", color= :grey45, #"#0087d7",
    font=:bold,
        )
        );
        tellheight=false,
        tellwidth=false,
        halign=0.99, 
        valign=1,
        fontsize=16,
        justification=:center, padding=(0,0,0,0)
        )

    Label(g_lay[1,1:3, Bottom()], rich(#"DEEP\nCUBE\n",
        rich("Explainable AI pipelines for big Copernicus data\nA Horizon 2020 research and innovation project",
        font=:regular,fontsize= 12*fs, color = 0.7*colorant"steelblue1"), font=:bold, color=colorant"dodgerblue3", fontsize= 16*fs,),
        tellwidth=false,tellheight=false,
        halign=0.0,
        #justification=:left,
        #valign=-0.5,
        padding=(0,0,0,70)
        )

    Label(fig[0,1:3], rich("Claire Robin, et. al. ",
        rich("Learning to forecast vegetation greenness at fine resolution over Africa with ConvLSTMs\n", color=0.5colorant"olivedrab2", font=:bold),
        rich("Artificial Intelligence for Humanitarian Assistance and Disaster Response.\nWorkshop at NeurIPS 2022, ",
        rich("https://doi.org/10.48550/arXiv.2210.13648", font=:regular), font=:bold, color=0.7*colorant"tan1"),
        color=:grey45, #1.35colorant"black",
        font=:bold, fontsize= 12*fs,
        ), tellwidth=false, tellheight=false,
        #valign=1.03,
        halign=0.75, justification=:left, padding=(0,0,10,10)
        )

    Box(fig[1,2], color=(0.85*colorant"tan1", 1), tellheight=false, tellwidth=false,
        halign=1.0, valign=0.93, width=130, height=30, cornerradius=10, strokevisible=false)

    Label(fig[1,2], rich("Target", color=:grey95, font=:bold, fontsize=16*fs);
        tellwidth=false, tellheight=false,
        halign=0.9, valign=0.94,
        justification=:right, padding=(0,0,10,10)
        )

    Box(fig[1,3], color=(0.5*colorant"steelblue1", 1), tellheight=false, tellwidth=false,
        halign=0.95, valign=0.93, width=130, height=30, cornerradius=10, strokevisible=false)

    Label(fig[1,3], rich("Prediction", color=:grey95, font=:bold, fontsize=16*fs);
        tellwidth=false, tellheight=false,
        halign=0.9, valign=0.94,
        justification=:right, padding=(0,0,10,10)
        )

    colgap!(fig.layout, 5)
    rowgap!(fig.layout, 2)
    rowsize!(fig.layout, 2, Auto(0.3))

    on(events(m.axis).mousebutton, priority = 2) do event
        if event.button == Mouse.left && event.action == Mouse.press
            # move marker
            plt, i = Makie.pick(m.axis)
            if plt == p
                idx_pos[] = i
                notify(idx_pos)
                #notify(name_var)
                return Consume(true)
            end
            #center!(ax.scene)
        end
        return Consume(false)
    end
    on(nan_menu.selection) do sel_color
        nan_color_o[] = sel_color
    end
    fig
end

with_theme(theme_light()) do
    #plot_cubes()
    #plot_cubes(; provider = CartoDB())
    plot_cubes(; provider = OpenTopoMap())
end