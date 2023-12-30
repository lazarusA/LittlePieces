# https://atmosphere.copernicus.eu/monitoring-ozone-layer
# https://ads.atmosphere.copernicus.eu/#!/search?text=ozone&keywords=((%20%22Product%20type:%20Reanalysis%22%20)%20AND%20(%20%22Spatial%20coverage:%20Global%22%20))
# https://codes.ecmwf.int/grib/param-db/?id=210203
# Ozone mass mixing ratio is the mass of ozone per kilogram of air.
using GLMakie
using YAXArrays, NetCDF, DimensionalData
using Healpix, NaNStatistics
using Statistics
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
function aggpixels!(m, data, indices)
    @sync begin
        for i ∈ eachindex(m)
           Threads.@spawn begin
                agg_val = nanmean(data[indices .== i])
                m.pixels[i] = agg_val
            end
        end
    end
end

function aggpixels(n, data, indices)
    return [nanmean(data[indices .== i]) for i ∈ 1:n]
end

function createPixMesh(resol, pointsperside; r = 1.0)
    n = resol.numOfPixels
    meshes = []
    points3d = []
    for pixidx in 1:n
        pts3d = boundariesRing(resol, pixidx, pointsperside, Float32)
        ptss = hcat(pts3d', pts3d[1,:])'
        pts3d = [r * Point3f(ptss[i,:]...) for i in axes(ptss, 1)]
        push!(points3d, pts3d)
        mf = triangle_mesh(points3d[1])
        mc = normal_mesh(pts3d, faces(mf))
        push!(meshes, mc)
    end
    return merge([meshes...])
end

# do a special marker
function createMarker(bottom_poly, f; h_b=1.0, h_t = 1.1)
    bottom_poly = h_b * bottom_poly
    top_poly = h_t .* bottom_poly
    top = GeometryBasics.Mesh(top_poly, f)
    bottom = GeometryBasics.Mesh(bottom_poly, f)
    combined = merge([top, bottom])
    nvertices = length(top.position)
    connection = Makie.band_connect(nvertices)
    m = normal_mesh(coordinates(combined), vcat(faces(combined), connection))
    return m
end

function createPixMesh3d(resol, pointsperside; h_b=1.0, h_t=1.05)
    n = resol.numOfPixels
    meshes = []
    points3d = []
    for pixidx in 1:n
        pts3d = boundariesRing(resol, pixidx, pointsperside, Float32)
        ptss = hcat(pts3d', pts3d[1,:])'
        pts3d = [Point3f(ptss[i,:]...) for i in axes(ptss, 1)]
        push!(points3d, pts3d)
        mf = triangle_mesh(points3d[1])
        mc = createMarker(pts3d, faces(mf); h_b, h_t);
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
function getPixIndicesNeighbours!(m, lon, lat, idx_dict, new_ind)
    for ϕ in lon, la in lat
        idx = ang2pix(m, θ(la), ϕ)
        idx_neig = getinterpolRing(m.resolution, θ(la), ϕ)[1]
        idx_dict[idx] = setdiff([idx_dict[idx]..., idx_neig...], [idx])
        push!(new_ind, idx)
    end
end

# Download data from here, snowfall 2023, December, all days, all hours.
# https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=form

ds = Cube(joinpath(@__DIR__, "ozone.nc"))
ds_ds = open_dataset(joinpath(@__DIR__, "ozone.nc"))

lon = Float64.(lookup(ds,:longitude)) * pi/180
lat = Float64.(lookup(ds,:latitude)) * pi/180
tempo = join.(split.(string.(lookup(ds, :Ti)), "T"), " ⋅ ")

all_surfaces = []
all_values= []
new_lon = range(0,2π, 1440)
new_lat = range(pi/2, -pi/2, 720)

for (i, r) in enumerate(range(1.005,1.1,25))
    # low res
    ds_data = replace(ds.data[:,:,1, 1], missing=>NaN)
    m, n, resol = createPixMap(5);
    indices = getPixIndices(m, lon, lat)
    aggpixels!(m, ds_data, indices);
    # high res
    new_m = [interpolate(m, θ(la), ϕ) for ϕ in new_lon, la in new_lat]
    m_new, n, resol = createPixMap(7);
    pix_mesh = createPixMesh(resol, 3; r=r);
    indices = getPixIndices(m_new, new_lon, new_lat)
    aggpixels!(m_new, new_m, indices);
    repeated_values = repeatValues(m_new, 3)
    push!(all_surfaces, pix_mesh)
    push!(all_values, repeated_values)
end

nan_mx = nanmaximum(m_new)
nan_mn = nanminimum(m_new)
sphere_w = uv_normal_mesh(Tesselation(Sphere(Point3f(0), 0.99), 128));

with_theme(theme_dark()) do
    fig = Figure(; size = (700,500),  backgroundcolor=0.5colorant"#232e41")
    ax = LScene(fig[1,1]; show_axis=false) 
    hm = mesh!(ax, all_surfaces[1];
        color=all_values[1],
        colorrange=(3.3589563f-6, 6.704774f-6),
        lowclip=:transparent,
        highclip=:yellow,
        colormap=resample_cmap(:lipari, 256, alpha=[x^2 for x in range(0,1,256)]),
        transparency=true,
        )
    for i in 2:15
        mesh!(ax, all_surfaces[i];
        color=all_values[i],
        colorrange=(3.3589563f-6, 6.704774f-6),
        lowclip=:transparent,
        highclip=:yellow,
        colormap=resample_cmap(:lipari, 256, alpha=[x^2 for x in range(0,1,256)]),
        transparency=true,
        )
    end

    mesh!(ax, sphere_w; color=img_r, shading=FastShading, transparency=false)
    Colorbar(fig[1,2], hm)
    fig
end








pix_mesh3d = createPixMesh3d(resol, 5; h_b=1.0, h_t=1.01);
repeated_values3d = repeatValues3d(m_new, 5);

img = load(joinpath(@__DIR__,"paper_world.png"));
img_r = circshift(img, (1,3600));
sphere_w = uv_normal_mesh(Tesselation(Sphere(Point3f(0), 1.0), 128));

with_theme(theme_minimal()) do
    fig = Figure(; size = (700,500), backgroundcolor=0.5colorant"#232e41"
        )
    ax = LScene(fig[1,1]; show_axis=false) 
    hm = mesh!(ax, pix_mesh3d;
        color=repeated_values3d,
        colorrange=(nan_mn, nan_mx),
        lowclip=:transparent,
        highclip=:yellow,
        colormap=resample_cmap(:Spectral, 256, alpha=[x^2 for x in range(0,1,256)]),
        shading=NoShading,
        transparency=true,
        )
    mesh!(ax, sphere_w; color=img_r, shading=FastShading, transparency=false)
    Colorbar(fig[1,2], hm)
    fig
end

# u_indices = unique(indices)
# idx_dict = Dict{Int, Array{Int}}()
# [idx_dict[u] = [] for u in eachindex(m)]; # surprise that this even works :D
# new_ind = []
# getPixIndicesNeighbours!(m, lon, lat, idx_dict, new_ind)
# pix2angRing(m.resolution, 1)
