using GLMakie
m = Rect3f(Vec3f(0), Vec3f(1))
function linkscenes!(connected_axes)
    scenes = [ax.scene for ax in connected_axes]
    locks = [Observable(false) for _ in scenes]
    for i in eachindex(scenes)
        cam = scenes[i].camera_controls
        onany(cam.eyeposition, cam.lookat, cam.upvector) do eyepos, lookat, up
            locks[i][] && return
            for j in eachindex(scenes)
                i == j && continue
                locks[j][] = true
                update_cam!(scenes[j], eyepos, lookat, up)
                locks[j][] = false
            end
        end
    end
end

fig = Figure()
ax1 = LScene(fig[1,1])
ax2 = LScene(fig[1,2])
mesh!(ax1, m; color = :red)
mesh!(ax2, m; color = :blue)
fig
connected_axes = [ax1, ax2]
linkscenes!(connected_axes)
fig