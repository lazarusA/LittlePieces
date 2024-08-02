using Healpix, NaNStatistics
using YAXArrays, NetCDF, DimensionalData
using GLMakie
using GLMakie.GeometryBasics
using GLMakie.FileIO
using Random

function createPixMap(order)
    resol = Resolution(2^order)
    m = HealpixMap{Float32, RingOrder}(resol.nside)
    n = m.resolution.numOfPixels
    return m, n, resol
end

θ(δ) = oftype(δ, π/2) - δ # co-latitude, δ-> latitude

"""
    getPixIndices(m, lon, lat)

m: a HealpixMap
lon: an array of longitude values in radians
lat: an array of latitude values in radians
"""
function getPixIndices(m, lon, lat)
    return [ang2pix(m, θ(la), ϕ) for ϕ in lon, la in lat]
end

"""
aggpixels!(m, n, data, indices)
"""
function aggpixels!(m, n, data, indices)
    @sync begin
        for i ∈ 1:n
           Threads.@spawn begin
            m.pixels[i] = nanmean(data[indices .== i])
           end
        end
    end
end

function aggpixels(n, data, indices)
    return [nanmean(data[indices .== i]) for i ∈ 1:n]
end

function createPixMesh(resol, pointsperside)
    n = resol.numOfPixels
    meshes = []
    points3d = []
    for pixidx in 1:n
        pts3d = boundariesRing(resol, pixidx, pointsperside, Float32)
        ptss = hcat(pts3d', pts3d[1,:])'
        pts3d = [Point3f(ptss[i,:]...) for i in axes(ptss, 1)]
        push!(points3d, pts3d)
        mf = triangle_mesh(points3d[1])
        mc = normal_mesh(pts3d, faces(mf))
        push!(meshes, mc)
    end
    return merge([meshes...])
end

# do a special marker
function createMarker(bottom_poly, f; h=1.0)
    top_poly = h .* bottom_poly
    top = GeometryBasics.Mesh(top_poly, f)
    bottom = GeometryBasics.Mesh(bottom_poly, f)
    combined = merge([top, bottom])
    nvertices = length(top.position)
    connection = Makie.band_connect(nvertices)
    m = normal_mesh(coordinates(combined), vcat(faces(combined), connection))
    return m
end

function createPixMesh3d(resol, m, pointsperside, α)
    n = resol.numOfPixels
    meshes = []
    points3d = []
    for pixidx in 1:n
        pts3d = boundariesRing(resol, pixidx, pointsperside, Float32)
        ptss = hcat(pts3d', pts3d[1,:])'
        pts3d = [Point3f(ptss[i,:]...) for i in axes(ptss, 1)]
        push!(points3d, pts3d)
        mf = triangle_mesh(points3d[1])
        mc = createMarker(pts3d, faces(mf); h = 1 + α * m[pixidx]);
        push!(meshes, mc)
    end
    return merge([meshes...])
end

function repeatValues(m, pointsperside)
    return repeat(m, inner=4*pointsperside + 1);
end

function repeatValues3d(m, pointsperside)
    return repeat(m, inner=4*pointsperside*2 + 2);
end

# Download data from here, snowfall 2023, December, all days, all hours.
# https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=form

ds = Cube(joinpath(@__DIR__, "snowfall.nc"))
ds_ds = open_dataset(joinpath(@__DIR__, "snowfall.nc"))
lon = Float64.(lookup(ds,:longitude)) * pi/180
lat = Float64.(lookup(ds,:latitude)) * pi/180
tempo = join.(split.(string.(lookup(ds, :Ti)), "T"), " ⋅ ")
ds_data = replace(ds.data[:,:,20], missing=>NaN)


m, n, resol = createPixMap(7);
pix_mesh = createPixMesh(resol, 3);
indices = getPixIndices(m, lon, lat);
aggpixels!(m, n, ds_data, indices);
repeated_values = repeatValues(m, 3)


img = load(joinpath(@__DIR__,"paper_world.png"));
img_r = circshift(img, (1,3600));
sphere_w = uv_normal_mesh(Tesselation(Sphere(Point3f(0), 0.99), 128));

with_theme(theme_minimal()) do
    fig = Figure(; size = (500,500),  backgroundcolor=0.5colorant"#232e41")
    ax = LScene(fig[1,1]; show_axis=false) 
    mesh!(ax, pix_mesh; color=repeated_values,
        colorrange=(0.000001,0.003),
        lowclip=:transparent,
        colormap=(:linear_kbc_5_95_c73_n256, 0.8),
        colorscale=log10
        )
    mesh!(ax, sphere_w; color=img_r, shading=FastShading, transparency=true)
    fig
end

m, n, resol = createPixMap(7);
ds_data = replace(ds.data[:,:, 1], missing=>NaN);
aggpixels!(m, n, ds_data, indices);

pix_mesh3d = Observable(createPixMesh3d(resol, m, 5, 30));
repeated_values3d = Observable(repeatValues3d(m, 5))
tempo_i = Observable(tempo[1])

set_theme!(theme_minimal())

fig = Figure(; size = (1280,720),  backgroundcolor=0.4colorant"#232e41") # 0.55colorant"#232e41" for png
axs = [LScene(fig[1,j]; show_axis=false) for j in 1:3]
hm = nothing
for j in 1:3
    hm = mesh!(axs[j], pix_mesh3d; color=repeated_values3d,
        colorrange=(0.000001,0.003),
        lowclip=:transparent,
        highclip=:orangered,
        colormap=Reverse(:Hiroshige),
        colorscale=log10
        )
end
[mesh!(axs[j], sphere_w; color=img_r, shading=FastShading, transparency=true) for j in 1:3]
Label(fig[0,1], rich("Global Snowfall: ", color=:orange, font=:bold, fontsize=20,
    rich("December 2023\n", color=:white),
    rich("Data: ERA5  -   ECMWF",
    color=:white, font=:regular, fontsize=14));
    tellwidth=false, halign=0.5, justification=:left
)
Label(fig[0,3], rich(rich("Visualization by ", color=:white), rich("Lazaro Alonso\n", color=:dodgerblue,
    rich("Created with ", color=:white), rich("Makie.jl", color=:orange)));
    tellwidth=false, halign=0.5,  justification=:left)
    
Box(fig[1,1:3], color=(0.25colorant"#232e41", 0.8), tellheight=false, tellwidth=false,
    valign=-0.015, halign=0.34, width=520, height=120, cornerradius=10)

Label(fig[1,1:3], rich("Motivated by the recent ", 
    rich("extreme snowfall event ", color=:dodgerblue), rich("in München Germany "),
    rich("where\nairplanes were frozen on the runway, "),
    rich("here we explore global snowfall patterns.\n"),
    rich("Since we are focusing mainly on the poles, ", color=:white),
    rich("the data presented here has been\n"),
    rich("aggregated into "), rich("equal-area regions using Healpix.jl. ", color=:orange, font=:bold),
    rich("This approach aims\nto provide a clearer interpretation "),
    rich("within these regions."), color=:white
    );
    tellheight=false, tellwidth=false,justification=:left,
    valign=0.01, halign=0.35
    )

Box(fig[1,1:3], color=(0.1colorant"#232e41", 0.95), tellheight=false, tellwidth=false,
    valign=0.97, halign=0.9, width=200, height=30, cornerradius=10)
Label(fig[1,1:3], @lift(rich($tempo_i, color=:white));
    tellheight=false, tellwidth=false,justification=:left,
    valign=0.96, halign=0.88)

cb=Colorbar(fig[1,2, Top()], hm, label=rich(" m of water equivalent", color=:white, font=:bold),
    vertical=false, ticklabelcolor=:white,
    tellwidth=false, tellheight=true, width=Relative(0.5), height=10,
    )
ax_snow = Axis(fig[0:1,1:3])
n_s = 55
Random.seed!(13)
plt1=scatter!(ax_snow, rand(n_s), rand(n_s); color=:white, marker= '❄', markersize=12)
plt2=scatter!(ax_snow, rand(n_s), rand(n_s); color=:white, marker= '❅', markersize=12)
plt3=scatter!(ax_snow, rand(n_s), rand(n_s); color=:white, marker= '❆', markersize=12)
limits!(ax_snow, 0.1,0.9,0.1,0.9)
hidedecorations!(ax_snow)
hidespines!(ax_snow)
#Makie.translate!(ax_snow.scene, 0, 0, -9900)
Makie.translate!(ax_snow.blockscene, 0, 0, -9900)
for (key, val) in ax_snow.interactions
    deregister_interaction!(ax_snow, key)
end
colgap!(fig.layout, 0)
rowgap!(fig.layout, 0)
fig

save(joinpath(@__DIR__, "../imgs/snowfall.png"), fig, update=false)

# the following is extremly slow!!!, because we are creating the meshes in each time step as well the aggregation needed.
record(fig, joinpath(@__DIR__,  "../imgs/snowfall.mp4"); framerate = 24, update=false) do io
    for i_time in eachindex(tempo)
        ds_data = replace(ds.data[:,:, i_time], missing=>NaN);
        aggpixels!(m, n, ds_data, indices);
        pix_mesh3d[] = createPixMesh3d(resol, m, 5, 30);
        repeated_values3d[] = repeatValues3d(m, 5)
        tempo_i[] = tempo[i_time]
        notify(pix_mesh3d)
        notify(repeated_values3d)
        notify(tempo_i)
        recordframe!(io)  # record a new frame
        @show i_time
    end
end