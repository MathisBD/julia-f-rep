
include(joinpath(dirname(Base.active_project()), "src", "csg", "csg.jl"))

using .Shapes, .Vec3s, .Nodes

# Create a Menger sponge of depth n>=0,
# with side length 1 and centered at the origin
function menger_sponge(n)
    @assert n >= 0
    if n == 0
        return Shapes.cube(Vec3(0.0, 0.0, 0.0), 1.0)
    else 
        m = menger_sponge(n-1)
        res = Shapes.empty
        for dx in -1.0:1.0
            for dy in -1.0:1.0
                for dz in -1.0:1.0
                    if abs(dx) + abs(dy) + abs(dz) > 1.0
                        res |= Shapes.translate(m, Vec3(dx, dy, dz))
                    end
                end
            end
        end
        return Shapes.scale(res, 1.0/3.0)
    end
end

# menger sponge 3 @ 0.1, 0.1, 0.1 ==> 315ms +- 64ms