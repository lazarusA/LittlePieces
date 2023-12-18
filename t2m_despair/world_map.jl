using CairoMakie
using CairoMakie.FileIO
using HDF5
CairoMakie.activate!()

f = h5open(joinpath(@__DIR__, "./world_xm.h5"), "r")
f110m = read(f, "world_110m")
lon = f110m["lon"]
lat = f110m["lat"]
pnts = Point2f.(lon, lat)

fig = Figure(figure_padding=0, size=(1440÷2,720÷2), backgroundcolor=:transparent)
ax = Axis(fig[1,1]; backgroundcolor=:transparent)
lines!(ax, pnts; color = :black, linewidth=0.025)
limits!(ax, -180,180,-90,90)
hidedecorations!(ax)
hidespines!(ax)
save(joinpath(@__DIR__, "../imgs/world_map.png"), fig)