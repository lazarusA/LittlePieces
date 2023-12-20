# code by Julius Krumbiegel and Lazaro Alonso
import ArchGDAL as AG
using GLMakie
# a similar output, different approach
# https://www.esri.com/arcgis-blog/products/arcgis-living-atlas/mapping/make-this-ai-inspired-topo-landscape-please/
# get data from https://dwtkns.com/srtm30m/, for this example look for Jena, Germany
dataset = AG.readraster(joinpath(@__DIR__, "N50E011.hgt"))
elev_band = AG.getband(dataset, 1)
elev_band_cut = rotr90(elev_band[1900:2550, 400:700]') # do a cut

cmap = ["#27755F","#5F9571","#E4DF60","#E2C14A","#BE591D"]
cmap = Reverse(:Hokusai1) #Reverse(:Hokusai1) #Reverse(:sandyterrain) #:starrynight #:sunset #Reverse(:sienna) # :Winter #:Bay #Reverse(:Homer2) #Reverse(:Hiroshige) #:thermal #:lipari
cmap = cgrad(resample_cmap(cmap, 100), scale=log2)
_, _, ct = contourf(elev_band_cut; levels=40,
    colormap = cmap)

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
fig = Figure(; figure_padding=0, size=(1200*fs,400*fs), backgroundcolor=0.5*cmap[1],
    fontsize= 10*fs)
ax = LScene(fig[1, 1], show_axis=false)
for (i, (poly, index)) in enumerate(zip(polys, level_indices))
    points = decompose(Point2f, poly.exterior)
    points_top = [Point3f(p..., index) for p in points]
    points_bottom = [Point3f(p..., 0) for p in points]

    faces = Makie.GeometryBasics.earcut_triangulate([points])
    msh = Makie.GeometryBasics.normal_mesh(
        Makie.to_ndim.(Point3f, points, index),
        faces
        )
    if !isempty(faces)
        mesh!(ax, msh; color = colors[i], shading=FastShading)
    end
    #lines!(ax, points_top; color=:black, linewidth=0.2, transparency=true) 
    band!(ax, points_bottom, points_top; color = colors[i], shading=FastShading)
end
Label(fig[0,1], rich(rich("Visualization by ", color=:white), rich("Lazaro Alonso & Julius Krumbiegel \n", color="white", font=:bold,
    rich("Created with ", color=:white), rich("Makie.jl", color=:greenyellow)));
    tellwidth=false, halign=0.01,  justification=:left, padding=(0,0,-10,10))

Label(fig[0,1], rich(rich("Topography of ", color=:white), rich("Jena, Germany \n", color="white", font=:bold,
    rich("Data ", color=:white), rich("30-Meter SRTM ", color=1.5colorant"orangered")));
    tellwidth=false, halign=0.99,  justification=:left, padding=(0,0,-10,10))
fig
cam = cameracontrols(ax.scene)
cam.lookat[] = Vec3f(277.91342, 137.40334, 9.034259)
cam.eyeposition[] = Float32[278.50653, 46.26349, 113.59202]
cam.upvector[] = Float32[-0.004905177, 0.7537961, 0.6570901]
cam.far[] = 2.2f0
cam.fov[] = 45.0f0
cam.near[] = 0.1f0
#display(fig; update=false)
save(joinpath(@__DIR__, "../imgs/elevation_jena.png"), current_figure(),
    update=false)