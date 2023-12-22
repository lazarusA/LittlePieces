using GLMakie
using JLD2
#GLMakie.activate!(float=true)

data_jld = jldopen(joinpath(@__DIR__, "sample_cubes.jld2"), "r")
sample_data = data_jld["target_prediction"]

function splat_values(data_array) # dims hard coded for this case
    return [data_array[i,j,k] for i in 1:10 for j in 1:128 for k in 1:128];
end
ndvi_points = [Point3f(i,j,k) for i in range(1,32,10) for j in 1:128 for k in 1:128];

fs=2
colormap = :Spectral
set_theme!()
fig = Figure(; size = (16*70*fs,9*50*fs), fontsize=12*fs, backgroundcolor=1.15colorant"gainsboro") #7colorant"#232e41"
axs = [LScene(fig[1,i]; show_axis=false) for i in 1:6]
axs2 = [LScene(fig[2,i]; show_axis=false) for i in 1:6]

[meshscatter!(axs[i], ndvi_points; marker=Rect3f(Vec3f(-0.5),Vec3f(1)),
        color= splat_values(sample_data[i][1]), colormap, shading=FastShading, colorrange=(0,1),
        lowclip=:grey9, nan_color= (:grey10, 0.05),
        markersize=Vec3f(0.8*3,0.99,0.99)) for i in 1:6]
[meshscatter!(axs2[i], ndvi_points; marker=Rect3f(Vec3f(-0.5),Vec3f(1)),
        color=splat_values(sample_data[i][2]), colormap, shading=FastShading, colorrange=(0,1),
        lowclip=:grey9, nan_color= (:grey10, 0.05),
        markersize=Vec3f(0.8*3,0.99,0.99)) for i in 1:6]
Label(fig[0,0:end], rich(rich("MAX PLANCK INSTITUTE\n", color=:black, font=:bold),
    rich("FOR BIOGEOCHEMISTRY", color=:grey8, fontsize=10*fs));
    tellwidth=false, halign=0.05,  justification=:right, padding=(0,0,-10,10)
    )

Label(fig[0,1:end], rich(rich("Visualization by ", color=:black), rich("Lazaro Alonso\n", color="black", font=:bold,
    rich("Created with ", color=:black), rich("Makie.jl", color=:orangered)));
    tellwidth=false, halign=0.99,  justification=:center, padding=(0,0,-10,10)
    )

Label(fig[3,0:end], rich("Claire Robin, Christian Requena-Mesa, Vitus Benson, Lazaro Alonso, Jeran Poehls, Nuno Carvalhais, Markus Reichstein\n",
    rich("Learning to forecast vegetation greenness at fine resolution over Africa with ConvLSTMs\n", color=0.5colorant"olivedrab2", font=:bold),
    rich("Artificial Intelligence for Humanitarian Assistance and Disaster Response.\nWorkshop at NeurIPS 2022, ",
    rich("https://doi.org/10.48550/arXiv.2210.13648", font=:regular), font=:regular, color=colorant"tan1"),
    color=1.35colorant"black", font=:bold, fontsize= 10*fs,
    ), tellwidth=false, halign=0.1,  justification=:left, padding=(0,0,-10,10)
    )

Label(fig[3,0:end], rich("DEEP\nCUBE\n", rich("Explainable AI pipelines for big Copernicus data\nA Horizon 2020 research and innovation project",
    font=:regular,fontsize= 10*fs, color = "steelblue1"), font=:bold, color=colorant"dodgerblue3", fontsize= 10*fs,),
    tellwidth=false, halign=1,  justification=:left, padding=(0,0,-10,10))

Colorbar(fig[1:2, 0],colormap=colormap, colorrange=(0,1), lowclip=:grey9,
        vertical=true, flipaxis=true, label="NDVI", labelrotation=0,
        #labelpadding=-20,
        height=Relative(0.3))
Label(fig[1,0], rich("Target", font=:bold, color=:tan1),  tellheight = false,tellwidth = true,)
Label(fig[2,0], rich("Prediction", font=:bold, color=:green),  tellheight = false,tellwidth = true,)
rowgap!(fig.layout,0)
colgap!(fig.layout,0)
fig
save(joinpath(@__DIR__, "../imgs/ndvi_samples.png"), fig)