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
fig = Figure(size=(1280,720), backgroundcolor="#2a1f12")
ax_bar = Axis(fig[2,1:3]; yaxisposition = :right,)
ax_img = Axis(fig[2,1:3]; aspect= DataAspect(), titlecolor=:silver, titlesize=12*fs,
    width= 150*fs, height= 120*fs, halign=0.0, valign=-0.2, tellwidth, tellheight)
axs = [LScene(fig[1,i], show_axis=false) for i in 1:3]

[surface!(axs[i], x,y,z; color = ds_data_mask, highclip=:black,
    lowclip=:grey8, colorscale = sc,
    colorrange = (-1,3), colormap, nan_color=:grey80) for i in 1:3]

[zoom!(axs[i].scene, cameracontrols(axs[i].scene), 0.65) for i in 1:3]
image!(ax_img, rotr90(img[10:end-10, 10:end-10]))

hm_o = barplot!(ax_bar, ds_data_series; color = ds_data_series,
    colorrange = (-1,2),  colormap=colormap_alpha, colorscale = sc, highclip=:black,
    lowclip=:grey8)
hm = barplot!(ax_bar, bar_lines; color = bar_lines,
    colorrange = (-1,2),  colormap, colorscale = sc,highclip=:black,
    lowclip=:grey8)
ax_bar.xticks = (slice_dates[2:end], tempo_s[slice_dates][2:end])
text!(ax_bar,Point2f(550,1.5), text=rich("Global mean", color = :dodgerblue))
translate!(hm_o, 0,0,100)
translate!(hm, 0,0,110)

xlims!(ax_bar, 1, length(ds_data_series))
ylims!(ax_bar, -1, 2.5)

rowsize!(fig.layout,2,Auto(0.2))
hidedecorations!(ax_img)
hidespines!(ax_img)
hideydecorations!(ax_bar)
hidexdecorations!(ax_bar, ticklabels=false)
ax_bar.xticklabelalign = (:right, :center)
# sl = Slider(fig[3, 1:3], range = 1:length(ds_data_series), startvalue = 1,
#     horizontal = true)
# connect!(t_step, sl.value)

cb=Colorbar(fig[2,3, Top()], hm, label=rich(" Data: ERA5 - ECMWF", color=:white),
    ticks = ([-1,0,1,2], ["-1 ᵒC", "0 ᵒC", "1 ᵒC", "2 ᵒC"]), vertical=false, ticklabelcolor=:white,
    tellwidth=false, tellheight=false, width=Relative(0.5), height=10,
    )

Label(fig[1:2,1:2], rich("Much like the scene captured in Géricault's ",
    rich(""""The Raft of the Medusa",\n""", color=:orangered, font=:bold),
    rich("our planet appears to be in chaos and besieged "),
    rich("by extraordinary events.\n", color=:white),rich("There's a pervasive sense of"),
    rich(" despair ", color=:white), rich("and "),
    rich("an urgent"), rich(" struggle for survival,\n", color=:white, font=:bold),
    rich("one that extends beyond our time, "),
    rich("echoing concerns for the well-being of\n"), rich("future generations.",
        color=:white, font=:bold)
    );
    tellheight, tellwidth,justification=:left, valign=0.13, halign=0.65
    )
Label(fig[0,1], rich("Global temperature anomalies\n", color=:orange, font=:bold, fontsize=20,
    rich("Reported as monthly anomalies relative to the 1950-1980 averages. ",
    color=:white, font=:regular, fontsize=14));
    tellwidth, halign=0, justification=:left
    )

Label(fig[0,3], rich(rich("Visualization by ", color=:white), rich("Lazaro Alonso\n", color=:dodgerblue,
    rich("Created with ", color=:white), rich("Makie.jl", color=:orange)));
    tellwidth, halign=0.5,  justification=:left)
Label(fig[0,2], rich("Don't look!\n", color =:dodgerblue,
    rich("It's worst than you think.", color=:orange));
    tellwidth, halign=0.5,  justification=:center)
colgap!(fig.layout,0)
rowgap!(fig.layout,0)
display(fig; update=false)

save(joinpath(@__DIR__, "../imgs/t2m_despair_w_x.png"), current_figure(),
    update=false)
# for i in [collect(1:39:876)..., 876]
#     t_step[] = i
#     save(joinpath(@__DIR__, "../imgs/t2m_despair_w_$(i).png"), current_figure(),
#         update=false)
# end

record(fig, joinpath(@__DIR__,  "../imgs/t2m_despair_z.mp4"); framerate = 32, update=false) do io
    for (idx,i) in enumerate(range(-1.5,4.5, 876))
        t_step[] = idx
        rotate!(axs[1].scene, i)
        rotate!(axs[2].scene, i)
        rotate!(axs[3].scene, i)
        recordframe!(io)  # record a new frame
    end
end
