module KdRaytracer
using Images, ..Vec3s, ..KdVoxels, ..Cameras, ..Tapes, ..TriDuals

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
    ishit = t_enter <= t_leave && 0.0 <= t_leave
    return (ishit, t_enter, t_leave) 
end

function raytrace_once(voxels :: KdVoxelGrid, ray :: Ray, t_curr :: Float64)
    level = 1
    node = voxels.node
    node_pos = voxels.world_pos
    cell_size = voxels.world_size / voxels.dims[level]

    # The time it takes to move one cell along a given axis.
    t_step = abs(cell_size / ray.dir)
    t_eps = 0.5 * minimum(t_step)
    @assert t_eps > 0.

    # The coordinates of the cell we are in.
    # Don't forget Julia indexing starts at 1.
    norm_pos = 1. + (ray(t_curr + t_eps) - node_pos) / cell_size
    coords = fmap(x -> floor(Int, x), norm_pos)

    # The time until we enter a new cell along a given axis.
    t_cross = t_curr + t_eps + 
        cell_size * (fmap(Float64, coords) + fmap(Float64, ray.dir >= 0.) - norm_pos) / ray.dir
    # A vector of +1 / -1, indicating how to change the coordinates when 
    # advancing along the ray.
    coords_step = Vec3s.ite(ray.dir >= 0., Vec3s.full(1), Vec3s.full(-1))

    # Step through the voxels one at a time along the ray.
    while 1 <= minimum(coords) && maximum(coords) <= voxels.dims[level]
        # We hit a voxel.
        if node.voxel_mask[coords.x, coords.y, coords.z] || (isa(node, InteriorNode) && node.node_mask[coords.x, coords.y, coords.z])
            return (true, t_curr)
        # Recurse in the child node.
        #elseif node <: InteriorNode && node.node_mask[coords.x, coords.y, coords.z]
        #    @assert level < length(voxels.dims)
        #    level += 1
        #    node = node.nodes[coords.x, coords.y, coords.z]
        #    node_pos += fmap(Float64, coords) * cell_size
        #    
        #    cell_size /= voxels.dims[level]
        #    t_step /= voxels.dims[level]
        #    t_eps /= voxels.dims[level]
        #    
        #    norm_pos = 1. + (ray(t_curr + t_eps) - node_pos) / cell_size
        #    coords = fmap(x -> floor(Int, x), norm_pos)
        #    @assert 1 <= minimum(coords) && maximum(coords) <= voxels.dims[level] 
        #
        #    t_cross = t_curr + t_eps + 
        #        cell_size * (fmap(Float64, coords) + fmap(Float64, ray.dir >= 0.) - norm_pos) / ray.dir
        # Step one cell forward.
        else
            # This is a vector of booleans with [true] where t_cross is minimal.
            # There can be several [true] coordinates.
            mask = t_cross == minimum(t_cross)
            @assert any(mask)

            t_curr = minimum(t_cross)
            t_cross += Vec3s.ite(mask, t_step, Vec3s.full(0.))
            coords += Vec3s.ite(mask, coords_step, Vec3s.full(0))
        end
    end

    # We didn't hit anything.
    return (false, t_curr)
end

# This uses the 'kd-restart' algorithm.
function raytrace(voxels :: KdVoxelGrid, ray :: Ray) :: Hit
    # Check the ray enters the voxel grid.
    (ishit, t_enter_grid, t_leave_grid) = raygrid_intersect(voxels, ray)
    if !ishit
        return Hit(-1.0)
    end

    # Raytrace.
    t_curr = max(t_enter_grid, 0.)
    fuel = 10000
    while t_curr < t_leave_grid - 0.01 && fuel > 0
        (ishit, t_next) = raytrace_once(voxels, ray, t_curr)
        @assert t_next >= t_curr
        if ishit
            return Hit(t_next)
        else 
            t_curr = t_next
        end

        fuel -= 1
    end

    return Hit(-1.)
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

function render!(img :: Array{RGB{N0f8}, 2}, camera :: Camera, voxels :: KdVoxelGrid, tape :: Tape)
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