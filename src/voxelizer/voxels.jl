module Voxels
using ..Vec3s

export VoxelGrid

struct VoxelGrid
    # The position / size of the grid in world coordinates.
    # The grid spans [world_pos] to [world_pos + world_size].
    world_pos :: Vec3{Float64} # world position of the low vertex of the grid
    world_cell :: Float64 # world size of a single voxel/cell of the grid
    world_size :: Float64 # world size of the entire grid
    # The number of voxels in each dimension.
    dim :: Int
    # The voxel grid has only one level (for the moment).
    mask :: BitArray{3}
end

function coords_to_world(grid :: VoxelGrid, coords :: Vec3{Int})
    return grid.world_pos + convert(Vec3{Float64}, coords) * grid.world_cell
end

function empty(world_pos :: Vec3{Float64}, world_size :: Float64, dim :: Int)
    world_cell = world_size / Float64(dim)
    mask = BitArray(fill(false, (dim, dim, dim)))
    return VoxelGrid(world_pos, world_cell, world_size, dim, mask)
end

end

