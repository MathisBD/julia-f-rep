module KdVoxelizer
using ..Vec3s, ..KdVoxels, ..Tapes, ..Intervals, ..Voxels

export voxelize, flatten


function voxelize_rec(tape :: Tape, node_pos :: Vec3{Float64}, cell_size :: Float64, dims :: Vector{Int}) :: Node
    d = dims[1]
    # Leaf Node.
    if length(dims) == 1
        node = empty_leaf_node(d)
        for z in 1:d
            for y in 1:d
                for x in 1:d
                    pos = node_pos + fmap(Float64, Vec3(x, y, z)) * cell_size 
                    res = Tapes.run(tape, pos.x, pos.y, pos.z)
                    node.voxel_mask[x, y, z] = res <= 0.
                end
            end
        end
        return node
    # Interior Node.
    else
        @assert length(dims) >= 2
        node = empty_interior_node(d)
        for z in 1:d
            for y in 1:d
                for x in 1:d
                    # First use interval evaluation.
                    pos = node_pos + fmap(Float64, Vec3(x, y, z)) * cell_size
                    res = Tapes.run(tape, 
                        Interval(pos.x, pos.x + cell_size), 
                        Interval(pos.y, pos.y + cell_size), 
                        Interval(pos.z, pos.z + cell_size))
                    # The voxel is full.
                    if res.high <= 0.
                        node.voxel_mask[x, y, z] = true
                    # The voxel is empty.
                    elseif res.low >= 0.
                        continue
                    # We need to create a child node.
                    else
                        node.node_mask[x, y, z] = true
                        node.nodes[x, y, z] = voxelize_rec(tape, pos, cell_size / d, dims[2:length(dims)])
                    end
                end
            end
        end
        return node
    end
end

function voxelize(tape :: Tape, world_pos :: Vec3{Float64}, world_size :: Float64, dims :: Vector{Int}) :: KdVoxelGrid
    @assert all(d -> d >= 1, dims)
    @assert length(dims) >= 1
    
    node = voxelize_rec(tape, world_pos, world_size / dims[1], dims)
    return KdVoxelGrid(world_pos, world_size, dims, node)
end

function flatten(voxels :: KdVoxelGrid) :: VoxelGrid

end

end