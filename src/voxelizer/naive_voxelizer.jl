module NaiveVoxelizer
using ..Vec3s, ..Voxels, ..Tapes

export voxelize!

function voxelize!(tape :: Tape, voxels :: VoxelGrid)
    for z in 1:voxels.dim
        for y in 1:voxels.dim
            for x in 1:voxels.dim
                pos = Voxels.coords_to_world(voxels, Vec3(x, y, z))
                res = Tapes.run(tape, pos.x, pos.y, pos.z)
                voxels.mask[x, y, z] = res <= 0.
            end
        end
    end
end

end