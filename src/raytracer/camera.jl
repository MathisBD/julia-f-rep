module Cameras
using ..Vec3s

export Camera, Ray

struct Camera
    pos :: Vec3{Float64}

    # (forward, up, right) should form an orthonormal basis.
    forward :: Vec3{Float64}
    up :: Vec3{Float64}
    right :: Vec3{Float64}
    
    # Horizontal field of view in radians
    fov_rad :: Float64
    # Aspect ratio (screen width / screen height) 
    aspect_ratio :: Float64
end

struct Ray
    origin :: Vec3{Float64}
    dir :: Vec3{Float64} # the direction should be normalized
end

@inline function (ray :: Ray)(t :: Float64) :: Vec3{Float64}
   return ray.origin + t * ray.dir 
end

# x and y are screen coordinates between 0.0 and 1.0
# (x-axis horizontal left-to-right and y-axis vertical top-to-bottom).
@inline function make_ray(camera :: Camera, x :: Float64, y :: Float64) :: Ray
    ux = tan(camera.fov_rad / 2) * camera.right
    uy = tan(camera.fov_rad / 2) / camera.aspect_ratio * camera.up
    dir = camera.forward + (2 * x - 1) * ux - (2 * y - 1) * uy
    return Ray(camera.pos, normalize(dir))
end

end