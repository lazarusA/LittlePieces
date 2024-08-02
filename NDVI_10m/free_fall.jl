using GLMakie
function the_fall(z, v, dt, g, restitution)
    z_new = z + v * dt - 0.5 * g * dt^2
    v_new = v - g * dt
    # Check for ground collision
    if z_new < 0
        # Bounce with damping
        z_new = -z_new * restitution
        v_new = -v_new * restitution
        # Reduce restitution for each bounce
        restitution *= 0.9
    end
    return z_new, v_new, restitution
end

function animate_fall!(zcor, duration, z, v, g)
    dt = 0.016  # Time step for updates (approx 60 FPS)
    t = 0.0
    restitution = fill(0.002, length(z))  # Coefficient of restitution for bounces
    while t < duration
        t += dt
        # Update position and velocity
        for (i, (zi, vi, gi, res_i)) in enumerate(zip(z, v, g, restitution))
            z_new, v_new, restitution_new = the_fall(zi, vi, dt, gi, res_i)
            z[i] = z_new
            v[i] = v_new
            restitution[i] = restitution_new
        end
        # Update the observable zcor
        zcor[] = z
        sleep(dt)  # Control the update frequency
        # Break the loop if the particle has essentially stopped
        if all(i -> i< 0.01, abs.(v)) && all(i -> i< 0.01, abs.(z)) # not a good line, it allocates!
            zcor[] .= 0  # Set to exactly 0 when very close
            break
        end
    end
end

function animate_fall_dt!(zcor, z, v, g, dt)
    restitution = fill(0.002, length(z))  # Coefficient of restitution for bounces
    # Update position and velocity
    for (i, (zi, vi, gi, res_i)) in enumerate(zip(z, v, g, restitution))
        z_new, v_new, restitution_new = the_fall(zi, vi, dt, gi, res_i)
        z[i] = z_new
        v[i] = v_new
        restitution[i] = restitution_new
    end
    # Update the observable zcor
    zcor[] = z
end


function update_z!(zcord, z)
    dt = 1e-6  # Time step for updates (approx 60 FPS)
    window_size = 711*8
    array_length = length(z)
    for i in 1:window_size:array_length
        end_index = min(i + window_size - 1, array_length)
        zcord.val[i:end_index] = z[i:end_index]
        zcord[] = zcord[]
        sleep(dt)  # Control the update frequency
    end
end

function update_z_window!(zcord, z)
    dt = 1e-6  # Time step for updates (approx 60 FPS)
    window_size = 711*8
    array_length = length(z)
    for i in 1:window_size:array_length
        end_index = min(i + window_size - 1, array_length)
        zcord.val[i:end_index] = z[i:end_index]
        zcord[] = zcord[]
    end
end

function generate_grid(xi, yi, xf, yf, n::Int, m::Int, x_gap, y_gap)
    width = xf - xi
    height = yf - yi
    
    dx = (width - (n-1) * y_gap) / n
    dy = (height - (m-1) * x_gap) / m  # Adjust dy to account for gaps
    
    grid = Rect3f[]
    
    for i in 0:n-1
        for j in 0:m-1
            x_min = xi + i * (dx + y_gap)
            y_min = yi + j * (dy + x_gap)  # Add gap to y_min calculation
            x_max = x_min + dx
            y_max = y_min + dy  # y_max is now y_min + dy (without gap)
            
            rect = Rect3f(Vec3f(x_min, y_min, 10.0f0), Vec3f(x_max - x_min, y_max - y_min, 20.0f0))
            push!(grid, rect)
        end
    end
    return grid
end

# # Initialize observable for the z-coordinate
# n = 10_000
# z_start = Observable(fill(10.0, n))
# xpos = 40*rand(n) .- 20
# ypos = 40*rand(n) .- 20

# ms = Observable(Vec3f(0.1, 0.1,0.25))
# # Create the figure and axis
# fig = Figure()
# ax = LScene(fig[1, 1])
# meshscatter!(ax, xpos, ypos, z_start; color = :grey25,
#     markersize = ms, marker = Rect3f(Vec3f(-0.5,-0.5, 0), Vec3f(1)))
# mesh!(ax, Rect3f(Vec3f(-40, -40, -0.1), Vec3f(80, 80, 0.1)))
# fig

# z_copy = z_start.val
# v_start = -2*rand(n)
# g_start = 3 .+ 10*rand(n)
# # Run the animation over 5 seconds with initial height 10, initial velocity 0, and gravity 9.8 m/s^2
# @time animate_fall!(z_start, 10, z_copy, v_start, g_start)
