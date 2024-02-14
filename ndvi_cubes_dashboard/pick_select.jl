using Random
using GLMakie
test_points = [Point2f(rand(2)...) for i in 1:10_000]
idx_pos = Observable(1)

mv_point = @lift(test_points[$idx_pos])
lbls = [randstring(3) for i in 1:10_000]

fig = Figure()
ax = Axis(fig[1,1])
p = scatter!(ax, test_points; markersize = 10, marker=:rect, color=:tan1,
    inspector_label = (self, i, p) -> lbls[i])
scatter!(ax, mv_point; color =(:white, 0.1), marker=:rect, markersize = 20, strokewidth=2,
    strokecolor=:red, inspector_label = (self, i, p) -> lbls[idx_pos[]]
    )
DataInspector(ax)
fig

on(events(ax).mousebutton, priority = 2) do event
    if event.button == Mouse.left && event.action == Mouse.press
        # move marker
        plt, i = Makie.pick(ax)
        if plt == p
            idx_pos[] = i
            notify(idx_pos)
            @show i
            return Consume(true)
        end
    end
    return Consume(false)
end
fig