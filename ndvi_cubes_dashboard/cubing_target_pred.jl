using GLMakie, YAXArrays
using Zarr, DimensionalData
using Colors
using TileProviders
using Tyler
using Extents
using MapTiles
using Statistics
using NaNStatistics

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
    yₘ = replace(yₘ, NaN=>missing)
    ŷₘ = replace(ŷₘ, NaN=>missing)

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

#ndvi_points_2 = [Point3f(i,j,k) for i in range(1,32,10), j in 1:128, k in 1:128];

# Observables
idx_pos = Observable(init_index)
mv_point = @lift(wm_points[$idx_pos])
# lift for cube name !
c_name = @lift(cubes_names[$idx_pos])

mn_cube = lift(c_name) do cube_n
    ndvi_data = c[ndvi=At("ndvi_target")]
    ndvi_data[cube_name=At(cube_n)].data[:,:,end:-1:1]
end

mn_cube_p = lift(c_name) do cube_n
    ndvi_data = c[ndvi=At("ndvi_pred")]
    ndvi_data[cube_name=At(cube_n)].data[:,:,end:-1:1]
end

μ_target = @lift([nanmean($mn_cube[:,:,k]) for k in 1:10])
μ_pred = @lift([nanmean($mn_cube_p[:,:,k]) for k in 1:10])

diff_t_p = @lift(abs.($μ_target .- $μ_pred))

ndvi_values = @lift([$mn_cube[i,j,k] for i in 1:10 for j in 1:128 for k in 1:128]);
ndvi_values_pred = @lift([$mn_cube_p[i,j,k] for i in 1:10 for j in 1:128 for k in 1:128]);

stats = @lift(getStats($mn_cube, $mn_cube_p))

nse = @lift($stats[:,:,1])
alpha = @lift($stats[:,:,2])
r_corr = @lift($stats[:,:,3])
beta = @lift($stats[:,:,4])


cmap = :Spectral
fs=1
ilon = minimum(lons) - 1
flon = maximum(lons) + 2
ilat = minimum(lats) - 1
flat = maximum(lats) + 0.2

extent = Extent(X = (ilon, flon), Y = (ilat, flat));

with_theme(theme_ggplot2()) do
    fig = Figure(; figure_padding=(15,15,5,0), size = (1400,900),
        backgroundcolor=1.15colorant"gainsboro",
        fontsize=16, fonts = (; regular="CMU Serif"))
    low_lay = GridLayout(fig[2,1:3])

    g_lay2 = GridLayout(low_lay[1,1:3])
    g_lay = GridLayout(low_lay[1,4:5])

    ax_inset = Axis(g_lay2[1,1:3], ylabel = "NDVI", xlabel="time step", #yaxisposition=:right
        )

    ax_map = Axis(fig[1,1])
    scatterlines!(ax_inset, μ_target, label = "spatial mean target",
        markersize=15,linestyle=:dash, color=:tan1
        )
    scatterlines!(ax_inset, μ_pred, label = "spatial mean prediction",
        marker=:rect, markersize=15, linestyle=:dash, color=0.7*colorant"steelblue1",
        )
    barplot!(ax_inset, 1:10, diff_t_p, width=0.2, color=:brown1, 
        label = "abs(target - prediction)")

    Legend(g_lay2[1,1:3], ax_inset, nbanks=3, tellwidth=false,tellheight = false,
        halign=0.5, valign=1.2)
    ax_inset.xticks = 1:10
    ylims!(ax_inset, 0, 0.6)
    xlims!(ax_inset, 0.5, 10.5)

    axs = [Axis(g_lay[1,i], aspect=1) for i in 1:4]

    hidedecorations!.(axs)
    hidespines!.(axs)
    hm1 = heatmap!(axs[1], nse, colorrange=(-0.5,1), colormap=:CMRmap, lowclip=:black)
    hm2 = heatmap!(axs[2], r_corr, colorrange=(-1,1), colormap=:Hiroshige,)
    hm3 = heatmap!(axs[3], alpha, colorrange=(0,3), colormap=:thermal, highclip=:gold)
    hm4 = heatmap!(axs[4], beta, colorrange=(-2.5,2.5), colormap=:Spectral_11,
        lowclip=:red4, highclip=:midnightblue)

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
        label=L"\alpha",
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

    m = Tyler.Map(extent; provider=CartoDB(), figure=fig,
        axis=ax_map, scale=3)

    p = scatter!(m.axis, wm_points; marker=:rect, markersize=10, color=(:tan1, 0.5),
        strokewidth=0.5, strokecolor=(:black,0.5))

    scatter!(m.axis, mv_point; marker=:rect, markersize=16, color=(:white, 0.1),
        strokewidth=2, strokecolor=:dodgerblue)

    hidedecorations!(ax_map)
    hidespines!(ax_map)

    ax = LScene(fig[1,2]; show_axis=false)

    plt = meshscatter!(ndvi_points; marker=Rect3f(Vec3f(-0.5),Vec3f(1)),
        color=ndvi_values, colormap=cmap, shading=FastShading, colorrange=(0,1),
        lowclip=:grey9, nan_color= (:grey10, 0.05),
        markersize=Vec3f(0.8*3,0.99,0.99))

    ax_p = LScene(fig[1,3]; show_axis=false)

    plt_p = meshscatter!(ax_p, ndvi_points; marker=Rect3f(Vec3f(-0.5),Vec3f(1)),
        color=ndvi_values_pred, colormap=cmap, shading=FastShading, colorrange=(0,1),
        lowclip=:grey9, nan_color= (:grey10, 0.05),
        markersize=Vec3f(0.8*3,0.99,0.99))

    Colorbar(fig[1,2, Bottom()], plt; vertical=false,
        label=@lift(rich("NDVI  ", rich("$($c_name)", font=:bold, color=:tan1))),
        width=Relative(0.3), tellwidth=false, halign=0)

    Label(fig[0,1:end], rich(rich("MAX PLANCK INSTITUTE\n", color=:black, font=:bold),
        rich("FOR BIOGEOCHEMISTRY", color=:grey8, fontsize=12*fs));
        tellwidth=false, halign=0.05,  justification=:right, padding=(0,0,10,10)
        )
    
    Label(g_lay[1,1:4, Bottom()], rich(rich("Visualization by ", color=:black), rich("Lazaro Alonso\n ", color="#0087d7", font=:bold,
        )
        );
        tellheight=false,
        tellwidth=false,
        halign=0.99, 
        valign=1,
        fontsize=16,
        justification=:center, padding=(0,0,10,10)
        )

    Label(fig[1,1:end, Top()], rich("DEEP\nCUBE\n",
        rich("Explainable AI pipelines for big Copernicus data\nA Horizon 2020 research and innovation project",
        font=:regular,fontsize= 12*fs, color = 0.7*colorant"steelblue1"), font=:bold, color=colorant"dodgerblue3", fontsize= 16*fs,),
        tellwidth=false,tellheight=false,
        halign=1,
        #justification=:left,
        #valign=-0.0001,
        padding=(0,0,20,0)
        )

    Label(fig[1,2:3], rich("Claire Robin, et. al. ",
        rich("Learning to forecast vegetation greenness at fine resolution over Africa with ConvLSTMs\n", color=0.5colorant"olivedrab2", font=:bold),
        rich("Artificial Intelligence for Humanitarian Assistance and Disaster Response.\nWorkshop at NeurIPS 2022, ",
        rich("https://doi.org/10.48550/arXiv.2210.13648", font=:regular), font=:bold, color=0.7*colorant"tan1"),
        color=1.35colorant"black", font=:bold, fontsize= 12*fs,
        ), tellwidth=false, tellheight=false,
        valign=1.1,
        halign=0.1, justification=:left, padding=(0,0,-10,10)
        )

    Box(fig[1,2], color=(0.85*colorant"tan1", 1), tellheight=false, tellwidth=false,
        halign=1.0, valign=0.93, width=130, height=30, cornerradius=10, strokevisible=false)

    Label(fig[1,2], rich("Target", color=:white, font=:bold, fontsize=16*fs);
        tellwidth=false, tellheight=false,
        halign=0.9, valign=0.94,
        justification=:right, padding=(0,0,10,10)
        )

    Box(fig[1,3], color=(0.5*colorant"steelblue1", 1), tellheight=false, tellwidth=false,
        halign=0.95, valign=0.93, width=130, height=30, cornerradius=10, strokevisible=false)

    Label(fig[1,3], rich("Prediction", color=:white, font=:bold, fontsize=16*fs);
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
                notify(name_var)
                return Consume(true)
            end
            #center!(ax.scene)
        end
        return Consume(false)
    end
    fig
end

#save(joinpath(@__DIR__, "../imgs/ndvi_dashboard_2x.png"), current_figure(), px_per_unit=2, update=false)
#save(joinpath(@__DIR__, "../imgs/ndvi_dashboard_2x_zoom.png"), current_figure(), px_per_unit=2, update=false)
#save(joinpath(@__DIR__, "../imgs/ndvi_dashboard_2x_zoom2.png"), current_figure(), px_per_unit=2, update=false)

# add TileProviders#main Tyler#master