module Voxels
using ..Vec3s

export VoxelGrid

struct VoxelGrid
    # The position / size of the grid in world coordinates.
    # The grid spans [world_pos] to [world_pos + world_size].
    world_pos :: Vec3{Float64} # world position of the low vertex of the grid
    world_cell :: Vec3{Float64} # world size of a single voxel/cell of the grid
    world_size :: Vec3{Float64} # world size of the entire grid
    # The number of voxels in each dimension.
    dim :: Vec3{Int}
    # The voxel grid has only one level (for the moment).
    voxels :: Array{Bool, 3}
end

function coords_to_world(grid :: VoxelGrid, coords :: Vec3{Int})
    return grid.world_pos + convert(Vec3{Float64}, coords) * grid.world_cell
end

function empty(world_pos :: Vec3{Float64}, world_size :: Vec3{Float64}, dim :: Vec3{Int})
    world_cell = world_size / convert(Vec3{Float64}, dim)
    voxels = fill(false, (dim.x, dim.y, dim.z))
    return VoxelGrid(world_pos, world_cell, world_size, dim, voxels)
end

end

