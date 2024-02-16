using GLMakie
time_scroll = Observable(0)

ndvi_points = [Point3f(i,j,k) for i in 1:10 for j in 1:3 for k in 1:3]
colors = Observable([i for i in 1:10 for j in 1:3 for k in 1:3])

scrolled_colors = @lift(circshift($colors, -9*($time_scroll*9)))

fig = Figure()
ax = LScene(fig[1,1])
meshscatter!(ax, ndvi_points; markersize = 0.5, color=scrolled_colors,
    colormap = :viridis,
    marker=Rect3f(Vec3f(-0.5), Vec3f(1)))
    
sg = SliderGrid(
    fig[1, 1, Bottom()],
    (label = "Time Scroller",range = 0:10, startvalue = 0, color_active=:black,
    color_inactive=(:grey9, 0.2), color_active_dimmed=:gold1,),
    width=Relative(0.5),
    tellheight = false)
connect!(time_scroll, sg.sliders[1].value)
fig