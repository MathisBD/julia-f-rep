module NaiveRaytracer
using Images, ..Vec3s, ..Voxels, ..Cameras, ..Tapes, ..TriDuals

export render!

const black = RGB{N0f8}(0, 0, 0)
const red = RGB{N0f8}(1, 0, 0)
const white = RGB{N0f8}(1, 1, 1)

# A negative [time] means there was no hit.
struct Hit 
    time :: Float64
end 

function raygrid_intersect(voxels, ray :: Ray)
    t_0 = (voxels.world_pos - ray.origin) / ray.dir
    t_1 = (voxels.world_pos + voxels.world_size - ray.origin) / ray.dir
    # The time at which the ray enters/leaves the grid along each axis.
    t_enter_axis = Vec3s.ite(ray.dir >= 0.0, t_0, t_1)
    t_leave_axis = Vec3s.ite(ray.dir >= 0.0, t_1, t_0)
    # The time at which the ray enters/leaves the grid globally.
    t_enter = max(t_enter_axis.x, t_enter_axis.y, t_enter_axis.z)
    t_leave = min(t_leave_axis.x, t_leave_axis.y, t_leave_axis.z)
    # Determine if the ray hit the grid
    hit = t_enter <= t_leave && 0.0 <= t_leave
    return (hit, t_enter, t_leave) 
end

function raytrace(voxels :: VoxelGrid, ray :: Ray) :: Hit
    # Check the ray enters the voxel grid.
    (hit, t_enter, _) = raygrid_intersect(voxels, ray)
    if !hit
        return Hit(-1.0)
    end
    
    # The time it takes to move one cell along a given axis.
    t_step = abs(voxels.world_cell / ray.dir)
    t_eps = 0.5 * minimum(t_step)
    @assert t_eps > 0.

    # The coordinates of the cell we are in.
    # Don't forget Julia indexing starts at 1.
    norm_pos = 1. + (ray(t_enter + t_eps) - voxels.world_pos) / voxels.world_cell
    coords = fmap(x -> floor(Int, x), norm_pos)
    
    # The time until we enter a new cell along a given axis.
    t_cross = t_enter + t_eps + 
        voxels.world_cell * (fmap(Float64, coords) + fmap(Float64, ray.dir >= 0.) - norm_pos) / ray.dir
    # A vector of +1 / -1, indicating how to change the coordinates when 
    # advancing along the ray.
    coords_step = Vec3s.ite(ray.dir >= 0., Vec3s.full(1), Vec3s.full(-1))

    # Step through the voxels one at a time along the ray.
    t = t_enter
    while 1 <= minimum(coords) && maximum(coords) <= voxels.dim
        # We hit something.
        if @inbounds voxels.mask[coords.x, coords.y, coords.z]
            return Hit(t)
        # Step one cell forward
        else
            # This is a vector of booleans with [true] where t_cross is minimal.
            # There can be several [true] coordinates.
            mask = t_cross == minimum(t_cross)
            @assert any(mask)

            t = minimum(t_cross)
            t_cross += Vec3s.ite(mask, t_step, Vec3s.full(0.))
            coords += Vec3s.ite(mask, coords_step, Vec3s.full(0))
        end
    end

    return Hit(-1.0)
end

function shade(tape :: Tape, ray :: Ray, hit :: Hit) :: RGB{N0f8}
    if hit.time < 0.
        return black
    else
        pos = ray(hit.time)
        x = TriDual(pos.x, 1., 0., 0.)
        y = TriDual(pos.y, 0., 1., 0.)
        z = TriDual(pos.z, 0., 0., 1.)
        output = Tapes.run(tape, x, y, z)
        
        grad = normalize(Vec3(output.dx, output.dy, output.dz))
        color = (grad + 1.0) / 2.0
        return RGB{N0f8}(color.x, color.y, color.z)
    end
end

function render!(img :: Array{RGB{N0f8}, 2}, camera :: Camera, voxels :: VoxelGrid, tape :: Tape)
    (screen_height, screen_width) = size(img)
    for screen_x in 1:screen_width
        for screen_y in 1:screen_height
            x = screen_x / Float64(screen_width)
            y = screen_y / Float64(screen_height)
            ray = Cameras.make_ray(camera, x, y)
            hit = raytrace(voxels, ray)
            color = shade(tape, ray, hit)
            img[screen_y, screen_x] = color
        end
    end
end

end