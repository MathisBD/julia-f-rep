module Cameras

struct Camera
    pos :: Vec3{Float64}

    # (forward, up, right) should form an orthonormal basis.
    forward :: Vec3{Float64}
    up :: Vec3{Float64}
    right :: Vec3{Float64}
    
    # Horizontal field of view in degrees
    fov_deg :: Float64
    # Aspect ratio (screen width / screen height) 
    aspect_ratio :: Float64
end

# screen_x and screen_y are coordinates between -1.0 and 1.0
# (x-axis horizontal and y-axis vertical).
# (-1, 1) ------ (1, 1)
#    |             |
#    |             |
# (-1,-1) ------ (1,-1)
function make_ray(screen_x :: Float64, screen_y :: Float64) :: Vec3{Float64}

end

end