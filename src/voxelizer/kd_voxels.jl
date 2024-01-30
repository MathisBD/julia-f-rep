module KdVoxels
using ..Vec3s

export Node, LeafNode, InteriorNode, KdVoxelGrid, empty_leaf_node, empty_interior_node

abstract type Node end

struct LeafNode <: Node
    voxel_mask :: BitArray{3}
end

struct InteriorNode <: Node
    voxel_mask :: BitArray{3}
    node_mask :: BitArray{3}
    nodes :: Array{Union{Nothing, Node}, 3}
end

function empty_leaf_node(dim :: Int)
    return LeafNode(BitArray(fill(false, (dim, dim, dim))))
end

function empty_interior_node(dim :: Int)
    return InteriorNode(
        BitArray(fill(false, (dim, dim, dim))),
        BitArray(fill(false, (dim, dim, dim))),
        fill(nothing, (dim, dim, dim)))
end

struct KdVoxelGrid
    # The position / size of the grid in world coordinates.
    # The grid spans [world_pos] to [world_pos + world_size].
    world_pos :: Vec3{Float64} # world position of the low vertex of the grid
    world_size :: Float64 # world size of the entire grid
    # The number of voxels in each dimension, for each level.
    dims :: Vector{Int}
    # The voxel grid has several levels.
    node :: Node
end

end

