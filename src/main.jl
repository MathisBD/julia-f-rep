include(joinpath(dirname(Base.active_project()), "src", "utils", "smooth_min_max.jl"))
include(joinpath(dirname(Base.active_project()), "src", "utils", "vec.jl"))
include(joinpath(dirname(Base.active_project()), "src", "utils", "tri_dual.jl"))
include(joinpath(dirname(Base.active_project()), "src", "utils", "interval.jl"))
include(joinpath(dirname(Base.active_project()), "src", "utils", "dot_graph.jl"))

include(joinpath(dirname(Base.active_project()), "src", "csg", "csg.jl"))
include(joinpath(dirname(Base.active_project()), "src", "voxelizer", "voxelizer.jl"))
include(joinpath(dirname(Base.active_project()), "src", "raytracer", "raytracer.jl"))

using ImageView, Images, Printf
using .Vec3s, .KdVoxels, .Cameras, .NaiveRaytracer, ..SmoothMinMax

# Create a Menger sponge of depth n>=0,
# with side length 1 and centered at the origin
function menger_sponge(n)
    @assert n >= 0 && n <= 4
    if n == 0
        return Shapes.cube(Vec3s.full(-0.5), 1.05)
    else 
        m = menger_sponge(n-1)
        res = Shapes.empty
        for dx in -1.:1.
            for dy in -1.:1.
                for dz in -1.:1.
                    if abs(dx) + abs(dy) + abs(dz) > 1.
                        res |= Shapes.translate(m, Vec3(dx, dy, dz))
                    end
                end
            end
        end
        return Shapes.scale(res, 1.01/3.0)
    end
end


# Create the shape.
shape = menger_sponge(1)
shape = Shapes.scale(shape, 12.)
shape = Shapes.rotateX(shape, pi / 4)
shape = Shapes.rotateY(shape, pi / 4)
shape = Shapes.rotateZ(shape, pi / 4)
shape &= -Shapes.sphere(Vec3(5.0, 0.0, 6.0), 7.0)

# Voxelize
tape = Tapes.node_to_tape(Nodes.constant_fold(shape))
println("Tape instruction count : $(length(tape.instructions))")

#dims = [4, 4, 4]
#voxels = Voxels.empty(Vec3s.full(-10.), 20., prod(dims))
#NaiveVoxelizer.voxelize!(voxels, tape)

#voxels = KdVoxelizer.voxelize(tape, Vec3s.full(-10.), 20., dims)
#voxels = KdVoxelizer.flatten(voxels)

# Create a target image.
width = 1200
height = 900
img = fill(RGB{N0f8}(0, 0, 0), (height, width))
# Create a camera looking down the z-axis.
camera = Camera(
    Vec3(0.0, 0.0, 30.0), # pos
    Vec3(0.0, 0.0, -1.0), # forward
    Vec3(0.0, 1.0, 0.0),  # up
    Vec3(1.0, 0.0, 0.0),  # right
    deg2rad(70.0),        # field-of-view in radians
    Float64(width) / Float64(height)) # aspect ratio

# Render the image
#render!(img, camera, voxels, tape)

# Display the image
#imshow(img)

function decomp(dim :: Int)
    @assert dim > 0

    MIN_FACTOR = 4
    MAX_FACTOR = 32

    if dim < MIN_FACTOR
        return []
    else
        results = []
        if dim <= MAX_FACTOR
            push!(results, [dim])
        end

        head = MIN_FACTOR
        while head <= dim && head <= MAX_FACTOR
            for tail in decomp(div(dim, head))
                push!(results, vcat([head], tail))
            end

            head *= 2
        end

        return results
    end
end

function benchmark_kd_voxelizer(dim :: Int)
    results = []
    dims_to_bench = decomp(dim)
    println("$(length(dims_to_bench)) configs to benchmark.")
    for dims in dims_to_bench
        println(">>> $dims")
        res = @benchmark KdVoxelizer.voxelize(tape, Vec3s.full(-10.), 20., $dims)
        push!(results, (dims, mean(res.times)))
    end

    sort!(results, by = x -> x[2])

    for (dims, time) in results
        @printf("%s\t ==> %.0fms\n", string(dims), time / 1e6)
    end
end

# Benchmarking kd-voxelizer (menger_sponge(1) with a sphere cut out) :
# [4, 8, 4, 4]     ==> 7002ms
# [8, 4, 4, 4]     ==> 7489ms
# [32, 4, 4]       ==> 7573ms
# [4, 4, 8, 4]     ==> 7927ms
# [16, 8, 4]       ==> 8090ms
# [8, 16, 4]       ==> 8866ms
# [4, 32, 4]       ==> 10466ms
# [4, 16, 8]       ==> 13371ms
# [16, 4, 8]       ==> 13945ms
# [8, 8, 8]        ==> 13989ms
# [4, 4, 4, 8]     ==> 14069ms
# [4, 8, 16]       ==> 25299ms
# [8, 4, 16]       ==> 26346ms
# [32, 16]         ==> 26916ms
# [4, 4, 32]       ==> 48322ms
# [16, 32]         ==> 48422ms