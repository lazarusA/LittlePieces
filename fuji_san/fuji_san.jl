using GLMakie, GLMakie.Colors, GLMakie.FileIO
using HDF5, Downloads
fid = h5open(joinpath(@__DIR__,"fuji.h5"), "r")
# Data subset from  https://hydro.iis.u-tokyo.ac.jp/~yamadai/MERIT_DEM/
topo = read(fid["fuji"])
link = "https://upload.wikimedia.org/wikipedia/commons/9/94/Umegawa_in_Sagami_province.jpg"
img = load(Downloads.download(link))

cmap = resample_cmap(:Hokusai3, 10)
cmap[9] = colorant"snow1"
cmap[10] = colorant"white"
tellwidth=false
tellheight=false
fs=3
with_theme(theme_dark()) do
    fig = Figure(; figure_padding=50, size= (600*fs,600*fs), backgroundcolor="#232e41")
    ax_o = Axis(fig[1,1]; title="『相州梅沢庄』⋅ 葛飾 北斎 ", titlecolor=:silver, titlesize=12*fs,
        width= 250*fs, height= 150*fs, halign=0.0, valign=0.915, tellwidth, tellheight)
    ax_t = Axis(fig[1,1])
    image!(ax_o, rotr90(img))
    c = 560
    for i in 1:560
        lines!(ax_t, topo[i, 50:end-30] .+ 15*c; linewidth = 0.85, overdraw = true,
            colorrange = (minimum(topo), maximum(topo)), color = topo[i, 50:end-30],
            colormap = cmap)
        c -= 1
    end
    limits!(ax_t, 1,289, 1,14000)
    Label(fig[1,1]; text=rich("富士山 ", color = 1.5colorant"firebrick", font=:bold,
        rich("(Fuji-san)", color = :grey85), fontsize=24*fs),
        tellwidth, tellheight, valign=1, halign=0.79)
    Label(fig[1,1]; text=rich("Created using 蒔絵 ", color = :white, font=:bold,
        rich(" [ Makie.jl ]", color = :gold), fontsize=12*fs),
        tellwidth, tellheight, valign=0.85, halign=0.8)
    Label(fig[1,1, Bottom()]; text=rich("Visualization by Lazaro Alonso", color=:white, fontsize=12*fs),
        tellwidth, halign=0.5)
    hidedecorations!.([ax_o, ax_t])
    hidespines!.([ax_o, ax_t])
    fig
end
save(joinpath(@__DIR__, "../imgs/fuji_san.png"), current_figure())
close(fid)
