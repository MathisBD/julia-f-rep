module NaiveRaytracer
using ..Vec3s, ..Voxels, ..Camera

export render!

function raytrace(voxels :: VoxelGrid, ray :: Vec3{Float64})
    return 0
end

function render!(camera :: Camera, grid :: VoxelGrid, img :: Array{Int, 2})
    (screen_width, screen_height) = size(img)
    for screen_x in 1:screen_width
        for screen_y in 1:screen_height
            x = 2.0 * (screen_x / screen_width) - 1.0
            y = 2.0 * (screen_y / screen_height) - 1.0
            ray = Camera.make_ray(camera, x, y)
            img[x, y] = raytrace(voxels, ray)
        end
    end
end

end