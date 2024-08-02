using GLMakie
using YAXArrays, NetCDF
using Dates
using DimensionalData
using NaNStatistics
using Statistics
using Distributions
using Random
using FlyThroughPaths

Random.seed!(1234)

include("free_fall.jl")

ds = open_dataset(joinpath(@__DIR__, "openEO_corr.nc"))
# ds = open_dataset(joinpath(@__DIR__, "openEO_milpas.nc"))

ndvi = ds["var"]
#ds["crs"]
xlon = lookup(ndvi, :x)
ylat = lookup(ndvi, :y)
xi, xmx = extrema(xlon)
yi, ymx = extrema(ylat)
xf = xmx - xi
yf = ymx - yi

boxGrid = generate_grid(xi, yi, xmx, ymx, 8, 4, 100, 50)

gpoints = [Point3f(i,j, 0) for i in xlon for j in ylat]

gpoints_4k = [Point3f(i,j, 4000) for i in xlon for j in ylat]

raw_data = ndvi[:,:,6].data

ndvi_r = readcubedata(ndvi)
axs = dims(ndvi_r)
data_all = ndvi_r.data[:,:,:] # read the data
ds = YAXArray(axs, data_all)

function mxYaxs(ds)
    g_ds = groupby(ds, :Ti => yearmonth)
    mx_g = maximum.(g_ds, dims=Ti)

    len_g = length(mx_g)
    window_size = 3
    last_start = len_g - window_size + 1

    ds_drei = map(1:window_size:last_start) do i
        ds_drei = cat(mx_g[i:i+window_size-1]..., dims=Ti)
        ds_drei_mx = mapslices(maximum, ds_drei, dims= "Ti");
        ds_drei_mx.data
    end
    return ds_drei
end
raw_data_all = mxYaxs(ds);
raw_data_all = cat(raw_data_all..., dims=3)

# gen labels
mabbr = monthabbr.(1:12)

# this should be in DD groupby
drei_months = ["Dec-Jan-Feb", "Mar-Apr-May", "Jun-Jul-Aug", "Sep-Oct-Nov"]
function generate_years_list(start_year::Int, end_year::Int)
    years = String[]
    push!(years, "$(start_year)-$(start_year+1)")
    for year in start_year+1:end_year
        push!(years, string.(repeat([year], 3))...)
        if year != end_year
            push!(years, "$(year)-$(year+1)")
        end
    end
    return years
end

years_str = generate_years_list(2015, 2024)
drei_months_str = repeat(drei_months, 9)

t_date = Observable(26)
zvals = @lift([raw_data_all[i,j,$t_date] for i in 1:length(xlon) for j in 1:length(ylat)])
#zvals = Observable(zvals)
# Define the points
b = Binomial(5000, 0.85)
gpoints_4k = [Point3f(i, j, rand(b)) for i in xlon for j in ylat]

xcord = [p[1] for p in gpoints_4k]
ycord = [p[2] for p in gpoints_4k]
zcord = Observable([p[3] for p in gpoints_4k])
cmap = resample_cmap(:fastie, 20)

# lights = [
#     AmbientLight(RGBf(0.75, 0.75, 0.75)),
#     #PointLight(RGBf(1, 1, 1), Point3f(0, 0, 1000), 5000),
#     # PointLight(RGBf(2, 0, 0), Point3f(xi, yi, 200), 1000),
#     # PointLight(RGBf(0, 2, 0), Point3f(-3,  3, 2), 1000),
#     # PointLight(RGBf(0, 0, 2), Point3f( 3,  3, 2), 1000),
#     # PointLight(RGBf(2, 2, 0), Point3f( 3, -3, 2), 1000),
# ]

function plotNDVI()
    cmap = resample_cmap(:fastie, 20)
    cmap2 = resample_cmap(Reverse(:Bay), 256)

    fig = Figure(;size = (1280,720), backgroundcolor=0.55colorant"#232e41"
        )
    ax = LScene(fig[1,1:3], show_axis=false,
        #scenekw = (; lights = lights)
        )
    meshscatter!(ax, gpoints;
        marker=Rect3f(Vec3f(-0.5, -0.5, 0), Vec3f(1)),
        markersize = Vec3f.(9.9, 9.9, 1),
        color = zvals,
        colormap = tuple.([:grey15, :grey50], 0.1),
        #colormap=cgrad(cmap[10:end], rev=false),
        colorrange = (0,1), lowclip=:grey15,
        transparency=true,
        )

    hm = meshscatter!(ax, xcord, ycord, zcord;
        marker=Rect3f(Vec3f(-0.5, -0.5, 0), Vec3f(1)),
        markersize = @lift(Vec3f.(9.9, 9.9, 150*($zvals))),
        color = zvals,
        colormap=cgrad(cmap[10:end], rev=false, alpha=1),
        colorrange = (0,1), lowclip=:grey15,
        specular=0.3,
        transparency=false)

    Label(fig[1,1], rich("From Space to Soil: ", color=1.25colorant"dodgerblue", font=:bold, fontsize=18,
        rich("NDVI", color=:orange, fontsize=16),
        rich("'s Role in Modern Farming\n", color=:white, fontsize=16),
        #rich("Maximum NDVI Dec-Jan-Feb 2023\n", color=:white, fontsize=16,),
        rich("Data: COPERNICUS - Sentinel-2_L2A",
        color=:white, font=:regular, fontsize=14));
        tellwidth=false,
        halign=0.5,
        valign=1.08,
        justification=:left,
        tellheight=false,
    )
    Label(fig[1,3], rich(rich("Visualization by ", color=:white),
        rich("Lazaro Alonso\n", color=1.25colorant"dodgerblue", font=:bold,
        rich("Created with ", color=:white), rich("Makie.jl", color=:orange)));
        tellwidth=false,
        halign=0.5,
        valign=1.08,
        tellheight=false,
        justification=:left)

    cb=Colorbar(fig[1,2], hm,
        label=@lift(
            rich("Maximum ", 
            rich("NDVI ", color = :orange), 
            rich("$(drei_months_str[$t_date]) ", color=1.25colorant"white" ),
            rich("$(years_str[$t_date])", color=:orange),
            color=:white, font=:bold)
            ),
        vertical=false, ticklabelcolor=:white,
        tellwidth=false, tellheight=false, width=Relative(0.95), height=10,
        valign=1, ticksize=10,
        )
    cb.ticks = 0:0.1:1

    Box(fig[1,1:3], color=(0.25colorant"#232e41", 0.8), tellheight=false,
        tellwidth=false, valign=-0.015, halign=1,
        width=650, height=140, cornerradius=10)

    Label(fig[1,1:3], rich("The Normalized Difference Vegetation Index ", 
        rich("(NDVI) ", color=:orange), rich("unveils the hidden vitality of "),
        rich("vegetation from space,\noffering a powerful "),
        rich("satellite-based tool  ", color=1.25colorant"dodgerblue", font=:bold),
        rich("with the potential to revolutionize agriculture. By analyzing how\nplants "),
        rich("reflect and absorb different wavelenghts of light, NDVI "),
        rich("delivers ", color=1.25colorant"dodgerblue", font=:bold),
        rich("crucial "),
        rich("insights ", color=1.25colorant"dodgerblue", font=:bold),
        rich("into "),
        rich("vegetation\nhealth and density. ", color=1.25colorant"dodgerblue", font=:bold),
        rich("The advent of ",),
        rich("new high-resolution 10m ", color=1.25colorant"dodgerblue", font=:bold),
        rich("data provides "),
        rich("unprecedented detail,\n",  color=1.25colorant"dodgerblue", font=:bold),
        rich("enabling farmers to precisely monitor field variations, ", color=:silver, font=:bold),
        rich("optimize resource use, and potentially\nboost yields, ",color=:silver, font=:bold),
        rich("contributing significantly to ",), 
        rich("global food security.", color=1.5*cmap2[1], font=:bold), color=:white
        );
        tellheight=false, tellwidth=false,justification=:left,
        valign=0.01, halign=1
        )
    cam = cameracontrols(ax.scene)
    cam.eyeposition[] = [228153.25, 2.4540862f6, 3986.6445]
    cam.lookat[] = [228166.67, 2.4592715f6, -1302.021]
    cam.upvector[] = [0.001849357, 0.71405077, 0.7000915]
    cam.fov[] = 45.0
    # zoom!(ax.scene, cam, 0.65)
    return fig, ax, cb
end
zcord = Observable([p[3] for p in gpoints_4k])

set_theme!(theme_ggplot2()) # set_theme!()
fig, ax, cb = plotNDVI()
display(fig; update=false)

z_copy = zcord.val
npts = length(z_copy)

v_start = -1000 .- rand(500:1500, npts).*abs.(randn(npts))
g_start = 10 .+ rand(60:100, npts).*abs.(randn(npts))
new_z = [raw_data_all[i,j,2] for i in 1:length(xlon) for j in 1:length(ylat)]
#animate_fall!(zcord, 15, z_copy, v_start, g_start)



view0 = capture_view(ax)
#@show view0
path = Path(view0)
path = path * ConstrainedMove(4,
    ViewState(
        eyeposition=[222491.6, 2.4566125f6, 2645.094],
        lookat=[228166.67, 2.4592715f6, -1302.021],
        upvector=[0.48257613, 0.22611324, 0.84616375],
        ), 
        :none, :constant
        )
path *= Pause(2)
path *= ConstrainedMove(4,
ViewState(
    eyeposition=[231409.55, 2.4583532f6, 1221.1691],
    lookat=[231362.4, 2.460477f6, -344.5451],
    upvector=[-0.013170958, 0.59313995, 0.8049916],
    ), 
    :none, :constant
    )
path *= Pause(2)
path *= ConstrainedMove(4,
    ViewState(
        eyeposition=[228652.38, 2.4549455f6, 1700.0756],
        lookat=[228556.58, 2.4594298f6, -1130.1426],
        upvector=[-0.011395099, 0.53351873, 0.84571147],
        ), 
        :none, :constant
        )

record(fig, joinpath(@__DIR__,  "../imgs/ndvi_from_space.mp4"), framerate = 24, update=false, px_per_unit=1.5) do io
    len_raw = last(axes(raw_data_all,3))
    for dt in 0:0.001:0.12
        animate_fall_dt!(zcord, z_copy, v_start, g_start, dt)
        recordframe!(io)  # record a new frame
    end
    for t in LinRange(0, 16, 150)
        set_view!(ax, path(t))
        recordframe!(io)  # record a new frame
    end

    pathPause = Path(capture_view(ax))
    pathPause *= ConstrainedMove(1, ViewState(), :none, :constant)
    pathPause *= Pause(1)

    pathEnd = Path(capture_view(ax))
    pathEnd *= ConstrainedMove(1, ViewState(
        eyeposition=[231554.45, 2.460527f6, 661.7063],
        lookat=[233090.66, 2.4604958f6, -448.77182],
        upvector=[0.58563685, -0.011960553, 0.8104853],),
        :none, :constant)

    window_size = 711*24
    array_length = length(new_z)

    for ts in axes(raw_data_all, 3)[27:end]
        new_z = [raw_data_all[i,j,ts] for i in 1:length(xlon) for j in 1:length(ylat)]
        window_size = ts == len_raw ? 711*8 : window_size
        for i in 1:window_size:array_length
            end_index = min(i + window_size - 1, array_length)
            zvals.val[i:end_index] = new_z[i:end_index]
            zvals[] = zvals[]
            recordframe!(io)  # record a new frame
        end
        t_date[] = ts

        for tpause in LinRange(0, 1, 50)
            set_view!(ax, pathPause(tpause))
            recordframe!(io)  # record a new frame
        end
        if ts==len_raw-1
            for tend in LinRange(0, 1, 50)
                set_view!(ax, pathEnd(tend))
                recordframe!(io)  # record a new frame
            end
            pathPause = Path(capture_view(ax))
            pathPause *= ConstrainedMove(1, ViewState(), :none, :constant)
            pathPause *= Pause(1)
        end
    end
    path_sum = Path(capture_view(ax))
    path_sum = path_sum * ConstrainedMove(1,
    ViewState(
        eyeposition=[228418.06, 2.452412f6, 6123.215],
        lookat=[228429.02, 2.4584748f6, -939.99384],
        upvector=[0.0013724946, 0.75879115, 0.65133256]
        ), 
        :none, :constant
        )
    path_sum *= Pause(1)

    for tpause in LinRange(0, 1, 150)
        set_view!(ax, path_sum(tpause))
        recordframe!(io)  # record a new frame
    end

    pathPause = Path(capture_view(ax))
    pathPause *= ConstrainedMove(0.01, ViewState(), :none, :constant)
    pathPause *= Pause(0.01)

    zvals[] .= NaN
    zvals[] = zvals[]

    cb.label[] = rich("Maximum ", rich("NDVI ", color = :orange), color=:white, font=:bold)
    t_date_new = Observable(1)

    Box(fig[1,1:3], color=:black, tellheight=false, tellwidth=false,
        valign=0.082, halign=0.20, width=200, height=40, cornerradius=10)
    Label(fig[1,1:3], @lift(
        rich("$(drei_months_str[$t_date_new]) ",
            rich("$(years_str[$t_date_new])", color=:orange),
            color=1.25colorant"white", font=:bold)
            ),
        tellheight=false, tellwidth=false,justification=:left,
        valign=0.1, halign=0.21)

    # do grid labels
    dx_r = 0.05:0.125:1 # range(0,1,10)
    dxx = dx_r[1]
    c = 1
    δd = 0.01
    for (k,y) in enumerate(2016:2023)
        if k <4
            dxx = dx_r[c] - δd/(3-(k-1))
        else
            dxx = dx_r[c] + δd/(8-(k-1))
        end
        Box(fig[1,1:3], color=(:black, 0.6), tellheight=false, tellwidth=false,
            valign=0.21, halign= dxx, width=75, height=20, cornerradius=10,
            )
        Label(fig[1,1:3], rich("$y", color=:orange, font=:bold),
            tellheight=false, tellwidth=false,justification=:left,
            valign=0.215, halign=dx_r[c])
        c +=1
    end

    for (k,valig) in enumerate([0.43, 0.6, 0.75, 0.9])
        Box(fig[1,1:3], color=:white, tellheight=false, tellwidth=false,
            valign=valig, halign= 1, width=85, height=20, cornerradius=10,
            )
        Label(fig[1,1:3], rich("$(drei_months[k])", color=:black, fontsize=12),
            tellheight=false, tellwidth=false,justification=:left,
            valign=valig, halign=0.995)
    end

    for k in 32:-1:1
        mesh!(ax, boxGrid[33-k]; color= Float32.(rotr90(raw_data_all[end:-1:1,:,33-k])),
            colormap=cgrad(cmap[10:end], rev=false, alpha=1),
            colorrange = (0,1), lowclip=:grey15, interpolate=false,)
        t_date_new[] = 33-k

        for tpause in LinRange(0, 1, 5)
            set_view!(ax, pathPause(tpause))
            recordframe!(io)  # record a new frame
        end
    end
end

save(joinpath(@__DIR__, "../imgs/ndvi_from_space.png"), fig, update=false, px_per_unit=1.5)

# mesh(boxGrid[1]; color= Float32.(raw_data_all[end:-1:1,:,1]),
#     colormap=cgrad(cmap[10:end], rev=false, alpha=1),
#     colorrange = (0,1), lowclip=:grey15, interpolate=false,)