using YAXArrays, Zarr, DimensionalData
using GLMakie.GeometryBasics
using GLMakie
using StatsBase
using GLMakie.FileIO
using Downloads
using Dates
GLMakie.activate!()

link = "https://upload.wikimedia.org/wikipedia/commons/thumb/1/15/JEAN_LOUIS_THÉODORE_GÉRICAULT_-_La_Balsa_de_la_Medusa_%28Museo_del_Louvre%2C_1818-19%29.jpg/640px-JEAN_LOUIS_THÉODORE_GÉRICAULT_-_La_Balsa_de_la_Medusa_%28Museo_del_Louvre%2C_1818-19%29.jpg"
img = load(Downloads.download(link));

# world map, coastlines
img_w = load(joinpath(@__DIR__, "../imgs/world_map.png"))
mask = Float32.(img_w .== img_w[1])
mask_z = mask .== 0
mask[mask_z] .= NaN;
mesh(Sphere(Point3f(0),1); color=mask, nan_color=:red)
# load data
# let's assume that you have a monthly dataset with lon,lat,time data, then doing the following will make sense :D
ds = Cube(joinpath(@__DIR__, "t2m_monthly.zarr/"))
ds_c = map(x-> x-273.15f0, ds)
t2m_ref = ds_c[Ti=Between(Date(1950), Date(1980))]

arr_ref = fill(NaN32, 1440, 720, 12)
for (i, m) in enumerate(["01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12"])
    m_inter = Date("1950-$(m)-01"):Year(1):Date("1980-01-01") |> collect
    t2m_ref = ds_c[Ti=At(m_inter)]
    t2m_μref = mapslices(mean, t2m_ref, dims="Time")
    arr_ref[:,:,i] .= t2m_μref.data[:,:]
end

#sphere = uv_normal_mesh(Tesselation(Sphere(Point3f(0), 1.5), 128))
tempo = lookup(ds_c, :Ti)
lentime = length(tempo)
slice_dates = range(1, lentime, step=lentime ÷ 7)
tempo_s = string.(Date.(tempo))

θ = range(0, π, 720)
φ = range(0, 2π, 1440)

t_step = Observable(876)
ds_data = lift(t_step) do t
    num_month = month(tempo[t])
    slice_t = ds_c[time=At(tempo[t])].data[:,:]' .- arr_ref[:,:,num_month]'
end

x = @lift([(1+$ds_data[i,j]/150)*cos(φ) * sin(θ) for (i,θ) in enumerate(θ), (j,φ) in enumerate(φ)])
y = @lift([(1+$ds_data[i,j]/150)*sin(φ) * sin(θ) for (i,θ) in enumerate(θ), (j,φ) in enumerate(φ)])
z = @lift([(1+$ds_data[i,j]/150)*cos(θ) for (i,θ) in enumerate(θ), (j,φ) in enumerate(φ)])

ds_data_mask = @lift($ds_data .* mask)

ds_data_series = Float32[]
for ti in eachindex(tempo)
    num_month = month(tempo[ti])
    anom = mean(ds_c[time=At(tempo[ti])].data[:,:] .- arr_ref[:,:,num_month])
    push!(ds_data_series, anom)
end

bar_lines = @lift(ds_data_series[1:$t_step])
fs=1.5
tellwidth, tellheight =false,false
colormap = :seaborn_icefire_gradient
n = 256
colormap = vcat(resample_cmap(:linear_kbc_5_95_c73_n256, n),
    resample_cmap(Reverse(:linear_kryw_5_100_c67_n256), n))

colormap_alpha = vcat(resample_cmap(:linear_kbc_5_95_c73_n256, n, alpha=0.1),
    resample_cmap(Reverse(:linear_kryw_5_100_c67_n256), n, alpha=0.1));
sc = ReversibleScale(x -> x > 0 ? x/2 : x, x -> x > 0 ? x*2 : x)

set_theme!(theme_dark())
fig = Figure(size=(1200,1200), backgroundcolor="#2a1f12")
ax = LScene(fig[1,1], show_axis=false)

surface!(ax, x,y,z; color = ds_data_mask, highclip=:black,
    lowclip=:grey8, colorscale = sc,
    colorrange = (-1,3), colormap, nan_color=:grey80)

#[zoom!(ax.scene, cameracontrols(ax.scene), 0.65) for i in 1:3]
display(fig; update=false)

# save(joinpath(@__DIR__, "../imgs/t2m_despair_obj.png"), current_figure(),
#     update=false)
