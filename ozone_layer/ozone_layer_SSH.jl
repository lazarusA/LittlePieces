using GLMakie
using GLMakie.FileIO
using GLMakie.GeometryBasics
using Statistics
using YAXArrays, NetCDF, DimensionalData
using Downloads
using SphereSurfaceHistogram
# ] add https://github.com/ffreyer/SphereSurfaceHistogram.jl.git

# download data from: 
# https://ads.atmosphere.copernicus.eu/cdsapp#!/dataset/cams-global-reanalysis-eac4?tab=form

# Inspiration and references:
# https://atmosphere.copernicus.eu/monitoring-ozone-layer
# https://www.unep.org/news-and-stories/story/rebuilding-ozone-layer-how-world-came-together-ultimate-repair-job
# https://www.youtube.com/watch?v=SDPxBZSXxpE

# Get data
p = @__DIR__
ds = Cube(joinpath(p, "ozone_07_12.nc"))

lon = Float64.(lookup(ds,:longitude)) * pi/180
lat = Float64.(lookup(ds,:latitude)) * pi/180
tempo = join.(split.(string.(lookup(ds, :Ti)), "T"), " ⋅ ")
ds_data = replace(ds.data[:,:,15, 1], missing=>NaN)

binners_sep_y = []
years = 2021:-1:2016
years_str = string.(years)
for year in years
    binners_sep =[]
    ds_c = Cube(joinpath(p, "ozone_09_$(year).nc"))
    # ds = mapslices(mean, ds_c, dims="level")[Ti=1]
    for l in 1:25
        binner = SSHAverager(10_000)
        ds_data = replace(ds_c.data[:,:,l, 1], missing=>NaN)
        append!(binner, pi/2 .- lat, lon, ds_data')
        push!(binners_sep, binner)
    end
    push!(binners_sep_y, binners_sep)
end

link = "https://eoimages.gsfc.nasa.gov/images/imagerecords/147000/147190/eo_base_2020_clean_720x360.jpg"
img = load(Downloads.download(link));
# img = load(joinpath(@__DIR__,"eo_base_2020_clean_720x360.jpg"));

img_r = circshift(img, (1,360));
sphere_w = uv_normal_mesh(Tesselation(Sphere(Point3f(0), 1), 64));
sphere_low = uv_normal_mesh(Tesselation(Sphere(Point3f(0), 1), 48));

# Set the goal number of bins
ddata = replace(ds.data[:,:,:, 1], missing=>NaN)

function get_binners(ddata, lat, lon; n_l=25, n_bins=10_000)
    binners = []
    for l in 1:n_l
        binner = SSHAverager(n_bins)
        ds_data = ddata[:,:,l]
        append!(binner, pi/2 .- lat, lon, ds_data')
        push!(binners, binner)
    end
    return binners
end

binners = get_binners(ddata, lat, lon)

n_l = 15
all_counts = [get_values(binner) for binner in binners[1:n_l]]
mn = minimum.(all_counts) |> minimum
mn = 1.75e-6
mx = maximum.(all_counts) |> maximum

# levels from 1hPa to 1000hPa
levels_hPa = ds.level
# barometric formula (simplest approach)
function h_altitude(P; R0=8.314462618, g=9.80665, M=0.02896968, T0=288.16, L=0.00976, P0=1013.25)
    return (T0/L) * (1 - exp(R0*L/(g*M)*(log(P/P0))))
end

h_alts = h_altitude.(levels_hPa)
rs = h_alts/maximum(h_alts)/15
tempo_i = Observable(tempo[1])

cmap = resample_cmap(Reverse(:Bay), 256);
colormap =  cgrad(cmap, alpha=[x^2.25 for x in range(0.01,0.75,256)]);
# generate visualization
meshes = []
for (r_i, binner) in enumerate(binners[1:n_l])
    m = vertex_mesh(binner, 1+rs[r_i])
    push!(meshes, m)
end
meshes = merge([meshes...])

m_values = vcat(get_values.(binners[1:n_l])...)
m_values_obs = Observable(m_values)

set_theme!(theme_ggplot2())

september=false

if september
    ddata = replace(ds.data[:,:,:, 497], missing=>NaN)
    binners = get_binners(ddata, lat, lon)
    m_values = vcat(get_values.(binners[1:n_l])...)
    m_values_obs = Observable(m_values)
    tempo_i = Observable(tempo[497])

end

fig = Figure(; size = (1280,720),  backgroundcolor=0.55colorant"#232e41") # 0.4colorant"#232e41" for mp4
ax = LScene(fig[1, 1], show_axis = !true)
ax_g = GridLayout(fig[1,2])
axs = []
for i in 2:3
    for j in 1:3
        push!(axs, LScene(ax_g[i,j]; show_axis=false))
    end
end

hm = mesh!(ax, meshes; color = m_values_obs, colorrange = (mn, mx), lowclip=:transparent,
    highclip=cmap[end], colormap, transparency = true)

mesh!(ax, sphere_w; color=img_r, shading=NoShading, transparency=false)

zoom!(ax.scene, cameracontrols(ax.scene), 0.55)
GLMakie.rotate!(ax.scene, Vec3f(-0.5, -1, 0.45), 3.0)

[mesh!(axs[k], sphere_low; color=img_r, shading=NoShading, transparency=false) for k in 1:6]

for k in 1:6
    binners_sep = binners_sep_y[k]
    hm = mesh!(axs[k], meshes; color = vcat(get_values.(binners_sep[1:n_l])...),
        colorrange = (mn, mx), colormap, lowclip=:transparent, highclip=cmap[end],
        transparency = true)
end

[zoom!(ax_i.scene, cameracontrols(ax_i.scene), 0.6) for ax_i in axs]
[GLMakie.rotate!(ax_i.scene, Vec3f(-0.5, -1, 0.45), 3.0) for ax_i in axs]
# label year
c = 1
for i in 2:3
    for j in 1:3
        Label(ax_g[i,j, Top()], years_str[c], tellwidth=false,
            tellheight=false, color=:white)
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
    rich("1970 levels.\n\n", color=1.25colorant"dodgerblue"),
    #rich("(see [1] and [2]).\n", color=:orange),
    rich("In this visualization we applied an equal area pixel aggregation over the\nsphere with "),
    rich("SphereSurfaceHistogram.jl", color=:lightseagreen),
    rich("."),
        color=:white
    );
    tellheight=false, tellwidth=false,justification=:left,
    valign=0.67, halign=0.5
    )

# Box(fig[1,1], color=(0.25colorant"#232e41", 0.95), tellheight=false, tellwidth=false,
#     valign=-0.015, halign=0.0, width=640, height=50, cornerradius=10)

# Label(fig[1,1], rich("[1] ",color=:orange, fontsize=12,
#     rich("https://www.unep.org/news-and-stories/story/rebuilding-ozone-layer-how-world-came-together-ultimate-repair-job\n", 
#         color=:white),
#     rich("[2] ", color=:orange),
#     rich("https://www.youtube.com/watch?v=SDPxBZSXxpE ", color=:white, fontsize=12)
#         );
#     tellheight=false, tellwidth=false,justification=:left,
#     valign=-0.0, halign=0.07
#     )
Label(fig[0,1:2], rich(rich("Visualization by ", color=:white),
    rich("Lazaro Alonso & Frederic Freyer\n ", color=1.25colorant"dodgerblue", font=:bold,
    rich("Created with ", color=:white), rich("Makie.jl", color=:lightseagreen)));
    tellwidth=false, halign=1.0,  justification=:left)
#Box(fig[1,1:2, Top()], color=:black, tellwidth=false, tellheight=true, width=Relative(0.35), height=10,)
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

if september
    save(joinpath(@__DIR__, "../imgs/ozone_hole.png"), fig, update=false)
end

record(fig, joinpath(@__DIR__,  "../imgs/ozone_laye_24_n.mp4"); framerate = 24, update=false) do io
    for i_time in eachindex(tempo) #data was corrupted after this time stamp
        ddata = replace(ds.data[:,:,:, i_time], missing=>NaN)
        new_binners = get_binners(ddata, lat, lon)
        m_new_values = vcat(get_values.(new_binners[1:n_l])...)
        m_values_obs[] = m_new_values
        tempo_i[] = tempo[i_time]
        notify(m_values_obs)
        notify(tempo_i)
        recordframe!(io)  # record a new frame
    end
end