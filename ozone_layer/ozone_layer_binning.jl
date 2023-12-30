################################################################################
### SphereSurfaceHistogram Extension
################################################################################

using GLMakie.GeometryBasics
using GLMakie.FileIO
using GLMakie
using Statistics
# may move this stuff into the package later...
using SphereSurfaceHistogram
# uses:
# partition_sphere1 or partition_sphere2
# partition_sphere_optim
# to_cartesian


struct SSHBinner2
    thetas::Vector{Float64}
    phi_divisions::Vector{Int64}

    # bin_many helpers
    phi_N_over_2pi::Vector{Float64}
    cumsum::Vector{Int64}
    theta_factor::Float64

    bins::Vector{Float64}
    count::Vector{Int64}
end

function SSHBinner2(N_bins::Int64; method = SphereSurfaceHistogram.partition_sphere2)
    thetas, phi_divisions = SphereSurfaceHistogram.partition_sphere_optim(4pi/N_bins, method=method)
    SSHBinner2(thetas, phi_divisions)
end

function SSHBinner2(thetas::Vector{Float64}, phi_divisions::Vector{Int64})
    # Precalculate N/2pi for bin_many!
    # use slightly more than 2pi to guarantee phi * phi_N_over_pi ∈ [0, 1)
    PI_PLUS = 2.0pi + 1e-14
    phi_N_over_2pi = phi_divisions ./ PI_PLUS

    # Also precalculate cumulative number of bins for each theta (from 0..pi)
    _cumsum = cumsum(phi_divisions)

    # to avoid conversions to Float64
    theta_factor = Float64(length(thetas)-1) / (pi + 1e-14)

    N_bins = sum(phi_divisions)
    bins = zeros(Float64, N_bins)
    count = zeros(Int64, N_bins)
    SSHBinner2(thetas, phi_divisions, phi_N_over_2pi, _cumsum, theta_factor, bins, count)
end

function Base.push!(B::SSHBinner2, longitude, latitude, value)
    theta_index = trunc(Int64, B.theta_factor * (pi/2 - latitude)) + 1
    phi_index = trunc(Int64, B.phi_N_over_2pi[theta_index] * longitude) + 1

    # Calculate position in bins
    l = theta_index > 1 ? B.cumsum[theta_index-1] : 0
    B.bins[l+phi_index] += value
    B.count[l+phi_index] += 1
    nothing
end

function bin_positions(B::SSHBinner2)
    points = Point3f[Point3f(0, 0, 1)]

    for i in 2:length(B.thetas)-2
        theta = 0.5(B.thetas[i] + B.thetas[i+1])
        N = B.phi_divisions[i]

        for k in 0:N-1
            phi = 2pi * (k + 0.5) / N
            push!(points, SphereSurfaceHistogram.to_cartesian(theta, phi))
        end
    end
    push!(points, Point3f(0, 0, -1))

    points
end

function generate_vertex_mesh(B::SSHBinner2; r = 1.0)
    faces = GLTriangleFace[]

    N = 0
    for i in 2:length(B.phi_divisions)
        top_idx = 1
        bot_idx = 1
        top = top_step = 1 / B.phi_divisions[i-1]
        bot = bot_step = 1 / B.phi_divisions[i]
        N2 = N + B.phi_divisions[i-1]

        while (bot < 0.99999) || (top < 0.99999)
            if top < bot
                top += top_step
                top_idx += 1
                i1 = N + mod1(top_idx, B.phi_divisions[i-1])
                i2 = N + top_idx - 1
                i3 = N2 + bot_idx
                push!(faces, GLTriangleFace(i1, i2, i3))
            else
                bot += bot_step
                bot_idx += 1
                i1 = N + top_idx
                i2 = N2 + bot_idx - 1
                i3 = N2 + mod1(bot_idx, B.phi_divisions[i])
                push!(faces, GLTriangleFace(i1, i2, i3))
            end
        end

        if i < length(B.phi_divisions)
            bot_idx += 1
            i1 = N + top_idx
            i2 = N2 + bot_idx - 1
            i3 = N2 + mod1(bot_idx, B.phi_divisions[i])
            push!(faces, GLTriangleFace(i1, i2, i3))
        end

        if i-1 > 1
            top_idx += 1
            i1 = N + mod1(top_idx, B.phi_divisions[i-1])
            i2 = N + top_idx - 1
            i3 = N2 + mod1(bot_idx, B.phi_divisions[i])
            top_idx = mod1(top_idx + 1, B.phi_divisions[i-1])
            push!(faces, GLTriangleFace(i1, i2, i3))
        end
        N = N2
    end

    # GeometryBasics.Mesh(GeometryBasics.meta(vertices; normals=ns), faces)
    normal_mesh(r * bin_positions(B), faces)
end

################################################################################
### Usage
################################################################################

# Grab data
using YAXArrays, NetCDF, DimensionalData
p = @__DIR__
ds = Cube(joinpath(p, "ozone.nc"))
ds_ds = open_dataset(joinpath(p, "ozone.nc"))

lon = Float64.(lookup(ds,:longitude)) * pi/180
lat = Float64.(lookup(ds,:latitude)) * pi/180
tempo = join.(split.(string.(lookup(ds, :Ti)), "T"), " ⋅ ")
ds_data = replace(ds.data[:,:,15, 1], missing=>NaN)

#GLMakie.volume(replace(ds.data[:,:,:, 1], missing=>NaN))

# Set the goal number of bins
binners = []
for l in 1:25
    binner = SSHBinner2(10_000)
    ds_data = replace(ds.data[:,:,l, 40], missing=>NaN)
    # bin
    for j in eachindex(lat)
        for i in eachindex(lon)
            push!(binner, lon[i], lat[j], ds_data[i, j])
        end
    end
    push!(binners, binner)
end

binners_sep_y = []
years = 2021:-1:2016
years_str = string.(years)
for year in years
    binners_sep =[]
    ds_c = Cube(joinpath(p, "ozone_09_$(year).nc"))
    # ds = mapslices(mean, ds_c, dims="level")[Ti=1]
    for l in 1:25
        binner = SSHBinner2(10_000)
        ds_data = replace(ds_c.data[:,:,l, 1], missing=>NaN)
        # bin
        for j in eachindex(lat)
            for i in eachindex(lon)
                push!(binner, lon[i], lat[j], ds_data[i, j])
            end
        end
        push!(binners_sep, binner)
    end
    push!(binners_sep_y, binners_sep)
end


img = load(joinpath(@__DIR__,"eo_base_2020_clean_720x360.jpg"));
# https://eoimages.gsfc.nasa.gov/images/imagerecords/147000/147190/eo_base_2020_clean_720x360.jpg
img_r = circshift(img, (1,360));
sphere_w = uv_normal_mesh(Tesselation(Sphere(Point3f(0), 1), 128));

sphere_low = uv_normal_mesh(Tesselation(Sphere(Point3f(0), 1), 48));

rs = range(1.1,1.005,25)

all_counts = [binner.bins ./ binner.count for binner in binners[1:15]]
mn = minimum.(all_counts) |> minimum
mx = maximum.(all_counts) |> maximum

# levels from 1hPa to 1000hPa
levels_hPa = ds_ds.level.val
# barometric formula (simples approach)
function h_altitude(P; R0=8.314462618, g=9.80665, M=0.02896968, T0=288.16, L=0.00976, P0=1013.25)
    return (T0/L) * (1 - exp(R0*L/(g*M)*(log(P/P0))))
end

h_alts = h_altitude.(levels_hPa)
rs = h_alts/maximum(h_alts)/15
tempo = join.(split.(string.(lookup(ds, :Ti)), "T"), " ⋅ ")
tempo_i = Observable(tempo[1])
cmap = resample_cmap(Reverse(:Bay), 256);
colormap =  cgrad(cmap, alpha=[x^2.25 for x in range(0.01,0.75,256)]);
# generate visualization

# https://eoimages.gsfc.nasa.gov/images/imagerecords/74000/74218/world.200412.3x5400x2700.png

with_theme(theme_ggplot2()) do
    fig = Figure(; size = (1280,720),  backgroundcolor=0.55colorant"#232e41") # 0.4colorant"#232e41" for mp4
    ax = LScene(fig[1, 1], show_axis = !true)
    ax_g = GridLayout(fig[1,2])
    #ax_text = Axis(ax_g[1,1:3])
    axs = []
    for i in 2:3
        for j in 1:3
            push!(axs, LScene(ax_g[i,j]; show_axis=false))
        end
    end

    hm=nothing
    for (r_i,binner) in enumerate(binners[1:15])
        m = generate_vertex_mesh(binner; r = 1+rs[r_i])
        hm = mesh!(ax, m; color = binner.bins ./ binner.count,
            colorrange = (mn - 0.1*mn, mx),
            lowclip=:transparent,
            highclip=cmap[end],
            colormap,
            transparency = true)
    end
    mesh!(ax, sphere_w; color=img_r, shading=NoShading, transparency=false)

    zoom!(ax.scene, cameracontrols(ax.scene), 0.625)
    GLMakie.rotate!(ax.scene, Vec3f(-0.5, -1, 0.45), 3.0)
    # for k in 1:6
    #     mesh!(axs[k], sphere_low; color=:white, shading=NoShading, transparency=false)
    # end
    for k in 1:6
        mesh!(axs[k], sphere_low; color=img_r, shading=NoShading, transparency=false)
        binners = binners_sep_y[k]
        for (r_i,binner) in enumerate(binners[1:15])
            m = generate_vertex_mesh(binner; r = 1+rs[r_i])
            hm = mesh!(axs[k], m; color = binner.bins ./ binner.count,
                colorrange = (mn - 0.1*mn, mx),
                lowclip=:transparent,
                highclip=cmap[end],
                colormap,
                transparency = true)
        end
    end
    [zoom!(ax_i.scene, cameracontrols(ax_i.scene), 0.6) for ax_i in axs]
    [GLMakie.rotate!(ax_i.scene, Vec3f(-0.5, -1, 0.45), 3.0) for ax_i in axs]
    c = 1
    for i in 2:3
        for j in 1:3
            Label(ax_g[i,j, Top()], years_str[c], tellwidth=false, tellheight=false, color=:white)
            c +=1
        end
    end

    Box(fig[1,1], color=:white, tellheight=false, tellwidth=false,
        valign=0.98, halign=0.025, width=200, height=30, cornerradius=10)
    Label(fig[1,1], @lift(rich($tempo_i, color=:black));
        tellheight=false, tellwidth=false,justification=:left,
        valign=0.97, halign=0.07)

    Box(fig[1,2], color=:black, tellheight=false, tellwidth=false,
        valign=0.63, halign=0.5, width=135, height=30, cornerradius=10)
    Label(fig[1,2], rich("09-01 ⋅ 00:00:00", color=:white);
        tellheight=false, tellwidth=false,justification=:left,
        valign=0.63, halign=0.5)

    Label(fig[0,1], rich("What happend to the ozone hole? Did it just go away?\n",
        rich("Data: CAMS global reanalysis (EAC4)  -   ECMWF",
        color=:white, font=:regular, fontsize=14), color=1.5*cmap[1], font=:bold, fontsize=16);
        tellwidth=false, halign=0.0, justification=:left
    )
    Label(fig[0,1:2], rich("""Ozone mass mixing ratio""", rich("""  kg ⋅ kg⁻¹ \n """, color=:white, fontsize=16),
        color=:lightseagreen, font=:bold, fontsize=20);
        tellwidth=false, halign=0.5, justification=:left
    )
    Box(ax_g[1,1:3], color=(0.25colorant"#232e41", 0.8), tellheight=false, tellwidth=false,
        valign=0.665, halign=0.34, width=500, height=170, cornerradius=10)

    Label(ax_g[1,1:3], rich("Ozone mass mixing ratio is the mass of ", 
        rich("ozone per kilogram of air. ", color=:lightseagreen),
        rich("So, what\nhappend to the ozone hole? "),
        rich("Well, "), rich("is still there, ", color=1.5cmap[1]), rich("however the trend over the \n"),
        rich("last five to ten years has been, ", color=:white),
        rich("NO INCREASE ", color=1.25colorant"dodgerblue", font=:bold),
        rich("[Montreal protocol] ", color = :orange),
        rich("in the\ninductive ozone hole. "),
        rich("Currently, estimates suggest that "), rich("it will be 2060 ", color=1.25colorant"dodgerblue", ),
        rich("before\nozone "), # https://www.youtube.com/watch?v=SDPxBZSXxpE
        rich("projections return ", color=1.25colorant"dodgerblue"),
        rich("to pre-"),
        rich("1970 levels ", color=1.25colorant"dodgerblue"),
        rich("(see [1] and [2]).\n", color=:orange),
        rich("In this visualization we applied an equal area pixel aggregation over the\nsphere with "),
        rich("SphereSurfaceHistogram.jl", color=:lightseagreen),
        rich("."),
            color=:white
        );
        tellheight=false, tellwidth=false,justification=:left,
        valign=0.67, halign=0.5
        )

    Box(fig[1,1], color=(0.25colorant"#232e41", 0.95), tellheight=false, tellwidth=false,
        valign=-0.015, halign=0.0, width=640, height=50, cornerradius=10)

    Label(fig[1,1], rich("[1] ",color=:orange, fontsize=12,
        rich("https://www.unep.org/news-and-stories/story/rebuilding-ozone-layer-how-world-came-together-ultimate-repair-job\n", 
            color=:white),
        rich("[2] ", color=:orange),
        rich("https://www.youtube.com/watch?v=SDPxBZSXxpE ", color=:white, fontsize=12)
            );
        tellheight=false, tellwidth=false,justification=:left,
        valign=-0.0, halign=0.07
        )
    Label(fig[0,1:2], rich(rich("Visualization by ", color=:white),
        rich("Lazaro Alonso & Frederic Freyer\n ", color=1.25colorant"dodgerblue", font=:bold,
        rich("Created with ", color=:white), rich("Makie.jl", color=:lightseagreen)));
        tellwidth=false, halign=1.0,  justification=:left)
    cb=Colorbar(fig[1,1:2, Top()],
        colorrange = (mn - 0.1*mn, mx),
        lowclip=:transparent,
        highclip=cmap[end],
        colormap = cgrad(resample_cmap(Reverse(:Bay), 256), 
            alpha=[15*x^3 for x in range(0.01,0.75,256)]),
        vertical=false, ticklabelcolor=:white,
        tellwidth=false, tellheight=true, width=Relative(0.35), height=10,
        )

    rowsize!(ax_g, 1, Auto(1.5))
    colgap!(ax_g, 0)
    rowgap!(ax_g, 0)
    colsize!(fig.layout, 1, Auto(1.5))
    display(fig, update=false)
end

