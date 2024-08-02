# using GLMakie
# cmap = resample_cmap(:fastie, 20)

# fig = Figure(size = (1280,600),)
# axs = [Axis(fig[i,j]) for j in 1:8 for i in 1:4]
# [axs[k].title = rich(years_str[k], color =:orange) for k in 1:4:32]
# [axs[k].ylabel = rich(drei_months[k], color =1.25*colorant"dodgerblue") for k in 1:4]

# [heatmap!(axs[k], raw_data_all[:,end:-1:1,k]; 
#     colormap=cgrad(cmap[10:end], rev=false, alpha=1),
#     colorrange = (0,1), lowclip=:grey15,) for k in 1:32]
# hidedecorations!.(axs, label=false)
# colgap!(fig.layout, 3)
# rowgap!(fig.layout, 3)
# fig




# fig = Figure()
# ax  = LScene(fig[1,1])
# meshscatter!(ax, some_points) # this one contains some good limits
# Tyler.Map3D(ax, my_region)
# fig

# cmap = resample_cmap(:fastie, 100)
# with_theme(theme_dark()) do
#     cmap = resample_cmap(:fastie, 20)
#     idx = Observable(1)
#     fig, ax, p = volume(@lift(ndvi[:,:, 1+ $idx:100 + $idx].data);
#         #colormap=:diverging_rainbow_bgymr_45_85_c67_n256,
#         #colormap=cmap[50:end],
#         colormap=cgrad(cmap[10:end], rev=false),
#         #colormap=:linear_protanopic_deuteranopic_kbw_5_98_c40_n256,
#         colorrange = (0,1), lowclip=:grey15)
#     Colorbar(fig[1,2], p)
#     sl = Slider(fig[end+1, 1:2], range = 1:525, startvalue = 1)
#     connect!(idx, sl.value)
#     fig
# end

# function create_intervals(xi, xf, n::Integer, gap)
#     if n < 1
#         return []
#     end
#     total_length = xf - xi
#     interval_length = (total_length - (n - 1) * gap) / n

#     intervals = []
#     current_start = xi
#     for _ in 1:n
#         interval_end = current_start + interval_length
#         push!(intervals, (current_start, interval_end))
#         current_start = interval_end + gap
#     end
#     return intervals
# end

# xr = create_intervals(first(xlon), last(xlon), 8, 5)
# yr = create_intervals(first(ylat), last(ylat), 4, 5)