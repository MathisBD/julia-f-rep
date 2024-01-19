module NaiveVoxelizer
using ..Vec3s, ..Voxels, ..Tapes

export voxelize!

function voxelize!(tape :: Tape, grid :: VoxelGrid)
    for x in 1:grid.dim.x
        for y in 1:grid.dim.y
            for z in 1:grid.dim.z
                pos = Voxels.coords_to_world(grid, Vec3(x, y, z))
                res = run(tape, pos.x, pos.y, pos.z)
                grid.voxels[x, y, z] = res <= 0
            end
        end
    end
end

end