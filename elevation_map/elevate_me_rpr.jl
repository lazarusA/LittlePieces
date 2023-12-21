# Currently this needs RPRMakie#master and is only tested on windows
using RPRMakie
using LinearAlgebra
import ArchGDAL as AG
using GLMakie

dataset = AG.readraster(joinpath(@__DIR__, "N50E011.hgt"))
elev_band = AG.getband(dataset, 1)
elev_band_cut = rotr90(elev_band[1900:2550, 400:700]') # do a cut

# helper funtions
gain(rgb::RGBf, gain) = RGBf(gain * rgb.r, gain * rgb.g, gain * rgb.b)
gain(rgba::RGBAf, gain) = RGBf(gain * rgba.r, gain * rgba.g, gain * rgba.b)

function azelrad(azim, elev, radius)
    x = radius * cosd(elev) * cosd(azim)
    y = radius * cosd(elev) * sind(azim)
    z = radius * sind(elev)
    Vec3f(x, y, z)
end

function kelvin_to_rgb(temp_kelvin)
    # Ensure temperature is within a valid range
    temp_kelvin = clamp(temp_kelvin, 1000, 40000)
    # Convert temperature to a scale of 0 to 1
    temp = temp_kelvin / 100
    # Calculate red
    red = temp <= 66 ? 255 : clamp(329.698727446 * ((temp - 60) ^ -0.1332047592), 0, 255)
    # Calculate green
    green = temp <= 66 ? clamp(99.4708025861 * log(temp) - 161.1195681661, 0, 255) : 
                         clamp(288.1221695283 * ((temp - 60) ^ -0.0755148492), 0, 255)
    # Calculate blue
    blue = temp >= 66 ? 255 : (temp <= 19 ? 0 : clamp(138.5177312231 * log(temp - 10) - 305.0447927307, 0, 255))
    # Convert to sRGB values (0 to 1 scale)
    return RGBf(red / 255, green / 255, blue / 255)
end

function arealight!(ax, matsys; azim, elev, radius, lookat, size, color)
    rect = Rect2f(-size/2, size)
    color_intensity_adjusted = lift(convert(Observable, color), convert(Observable, size)) do color, size
        Makie.to_color(color) ./ *(size...)
    end
    m = mesh!(ax, rect, material = RPR.EmissiveMaterial(matsys), color = color_intensity_adjusted)

    q = Makie.rotation_between(Vec3f(0, 0, 1), -azelrad(azim, elev, 1))
    Makie.rotate!(m, q)
    Makie.translate!(m, azelrad(azim, elev, radius) + lookat)
    return m
end

function saveRPR(filename::String, imageOut; new_size = (1200,400),  backgroundcolor=:grey8, fs=1)
    fig = Figure(; figure_padding=0, size= new_size,  backgroundcolor= backgroundcolor, fontsize=10*fs)
    ax = Axis(fig[1,1], aspect = DataAspect())
    image!(ax, rotr90(imageOut))
    Label(fig[0,1], rich(rich("MAX PLANCK INSTITUTE\n", color=:white, font=:bold),
    rich("FOR BIOGEOCHEMISTRY", color=:white, fontsize=7*fs));
        tellwidth=false, halign=0.01,  justification=:right, padding=(0,0,-10,10)
        )

    Label(fig[0,1], rich(rich("Visualization by ", color=:white), rich("Lazaro Alonso & Julius Krumbiegel \n", color="white", font=:bold,
    rich("Created with ", color=:white), rich("Makie.jl", color=:greenyellow)));
        tellwidth=false, halign=0.5,  justification=:center, padding=(0,0,-10,10)
        )

    Label(fig[0,1], rich(rich("Topography of ", color=:white), rich("Jena, Germany \n", color="white", font=:bold,
        rich("Data ", color=:white), rich("30-Meter SRTM ", color=1.5colorant"orangered")));
        tellwidth=false, halign=0.99,  justification=:left, padding=(0,0,-10,10))
    hidedecorations!(ax)
    hidespines!(ax)
    save(joinpath(@__DIR__, "../imgs/$(filename).png"), fig)
end

# RPRMakie.activate!(resource=RPR.RPR_CREATION_FLAGS_ENABLE_CPU,
#         plugin=RPR.Northstar, iterations=10)
    
RPRMakie.activate!(resource=RPR.RPR_CREATION_FLAGS_ENABLE_GPU0,
    plugin=RPR.Northstar, iterations=1500)

cmap = ["#27755F","#5F9571","#E4DF60","#E2C14A","#BE591D"]
cmap = Reverse(:Hokusai1) #Reverse(:Hokusai1) #Reverse(:sandyterrain) #:starrynight #:sunset #Reverse(:sienna) # :Winter #:Bay #Reverse(:Homer2) #Reverse(:Hiroshige) #:thermal #:lipari
cmap = cgrad(resample_cmap(cmap, 100), scale=log2)
_, _, ct = contourf(elev_band_cut; levels=40, colormap = cmap)
        
polys = ct.plots[1][1][]
level_values = ct.plots[1].color[]
level_indices = let 
    i = 0
    prev = nothing
    indices = Int[]
    for val in level_values
        if val != prev
            i += 1
        end
        push!(indices, i)
        prev = val
    end
    indices
end
colors = Makie.numbers_to_colors(ct.plots[1].color[], ct.plots[1])
fs = 2
lights = [
    # EnvironmentLight(1.0, [RGBf(0, 0, 0);;]),
    AmbientLight(kelvin_to_rgb(20_000)),
    DirectionalLight(gain(kelvin_to_rgb(3000), 2), -azelrad(270, 60, 1)),
]

sc = Scene(camera = cam3d!, size=(1200*fs,400*fs),  backgroundcolor=0.5*cmap[1])
sc.lights = lights

screen = RPRMakie.Screen(sc)
matsys = screen.matsys
context = screen.context
mat = RPR.DiffuseMaterial(matsys)

for (i, (poly, index)) in enumerate(zip(polys, level_indices))
    points = decompose(Point2f, poly.exterior)
    points_top = [Point3f(p..., (index)) for p in points]
    points_bottom = [Point3f(p..., 0) for p in points]

    faces = Makie.GeometryBasics.earcut_triangulate([points])

    msh = Makie.GeometryBasics.normal_mesh(
        Makie.to_ndim.(Point3f, points, index),
        faces
    )
    if !isempty(faces)
        mesh!(sc, msh, color = colors[i])
    end

    band!(sc, points_bottom, points_top, color = gain(colors[i], 1))
end
lims = Makie.data_limits(sc)

cam = cameracontrols(sc)
cam.lookat[] = Vec3f(277.91342, 137.40334, 9.034259)
cam.eyeposition[] = Float32[278.50653, 46.26349, 113.59202]
cam.upvector[] = Float32[-0.004905177, 0.7537961, 0.6570901]
cam.far[] = 1000f0
cam.fov[] = 45.0f0
cam.near[] = 0.1f0
GLMakie.activate!(float = true)

# sc

display(screen, sc)

##

# @time display(screen)

context, task = RPRMakie.replace_scene_rpr!(sc, screen)

##

# only display(screen, sc) honors camera settings plus fov for some reason

#save to file  with credits

display(screen, sc)
imageOut = colorbuffer(screen);
saveRPR("elevation_jena_rpr", imageOut; new_size=(1200*fs,430*fs), backgroundcolor = 0.5*cmap[1], fs=fs)
