using GLMakie
using Dates
using YAXArrays, Zarr

ds = open_dataset(joinpath(@__DIR__, "ndvi_clean_4534_cubes.zarr/"))
tempo = ds["tempo_target"]

t_xticks = tempo[cube_name=At("30RTU1877")].data[:]
t_str = string.(t_xticks)
d_str = string.(Date.(t_str, DateFormat("yyyymmdd")))
d_str_o = Observable(d_str)

fig, ax, obj = scatter(1:10)
ax.xticklabelrotation = π / 4
ax.xticklabelalign = (:right, :center)
fig
lift(d_str_o) do n_str
    ax.xticks = (1:10, n_str)
end
fig
