include(joinpath(dirname(Base.active_project()), "src", "utils", "smooth_min_max.jl"))
include(joinpath(dirname(Base.active_project()), "src", "utils", "vec.jl"))
include(joinpath(dirname(Base.active_project()), "src", "utils", "tri_dual.jl"))
include(joinpath(dirname(Base.active_project()), "src", "utils", "interval.jl"))
include(joinpath(dirname(Base.active_project()), "src", "utils", "dot_graph.jl"))

include(joinpath(dirname(Base.active_project()), "src", "csg", "csg.jl"))
include(joinpath(dirname(Base.active_project()), "src", "voxelizer", "voxelizer.jl"))
include(joinpath(dirname(Base.active_project()), "src", "raytracer", "raytracer.jl"))

using ImageView, Images, Printf, CUDA, StaticArrays
using .Vec3s, .KdVoxels, .Cameras, .NaiveRaytracer, .Tapes

# Create a Menger sponge of depth n>=0
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

function main_old()
    # Create the shape.
    shape = menger_sponge(2)
    shape = Shapes.scale(shape, 12.)
    shape = Shapes.rotateX(shape, pi / 4)
    shape = Shapes.rotateY(shape, pi / 4)
    shape = Shapes.rotateZ(shape, pi / 4)
    shape &= -Shapes.sphere(Vec3(5.0, 0.0, 6.0), 7.0)

    # Voxelize
    tape = node_to_tape(Nodes.constant_fold(shape))
    println("Tape instruction count : $(length(tape.instructions))")

    dims = [8, 4, 4, 4]
    #voxels = Voxels.empty(Vec3s.full(-10.), 20., prod(dims))
    #NaiveVoxelizer.voxelize!(voxels, tape)

    voxels = KdVoxelizer.voxelize(tape, Vec3s.full(-10.), 20., dims)
    voxels = KdVoxelizer.flatten(voxels)

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
    render!(img, camera, voxels, tape)

    # Display the image
    imshow(img)
end

const MAX_SLOT_COUNT = 32

function voxelize_kernel(
    dim :: Int,
    grid_pos :: Vec3{Float64},
    grid_size :: Float64,
    instructions_d :: CuDeviceArray{Instruction}, 
    constants_d :: CuDeviceArray{Float64},
    voxels_d :: CuDeviceArray{Bool, 3})

    x = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    y = (blockIdx().y - 1) * blockDim().y + threadIdx().y
    z = (blockIdx().z - 1) * blockDim().z + threadIdx().z
    #thread = (threadIdx().z - 1) * blockDim().y * blockDim.x + 
    #         (threadIdx().y - 1) * blockDim().x + 
    #         threadIdx().x 
    
    slots = zeros(MVector{MAX_SLOT_COUNT, Float64})
    slots[1] = (x - 1) * grid_size / Float64(dim) + grid_pos.x 
    slots[2] = (y - 1) * grid_size / Float64(dim) + grid_pos.y
    slots[3] = (z - 1) * grid_size / Float64(dim) + grid_pos.z
    
    @inbounds for inst in instructions_d
        if inst.op == Copy
            slots[inst.out_slot] = slots[inst.in_slotA]
        elseif inst.op == LoadConst
            slots[inst.out_slot] = constants_d[inst.in_slotA]
        elseif inst.op == Sin
            slots[inst.out_slot] = sin(slots[inst.in_slotA])
        elseif inst.op == Cos
            slots[inst.out_slot] = cos(slots[inst.in_slotA])
        elseif inst.op == Exp
            slots[inst.out_slot] = exp(slots[inst.in_slotA])
        elseif inst.op == Neg
            slots[inst.out_slot] = -slots[inst.in_slotA]
        elseif inst.op == Sqrt
            slots[inst.out_slot] = sqrt(slots[inst.in_slotA])
        elseif inst.op == Add
            slots[inst.out_slot] = slots[inst.in_slotA] + slots[inst.in_slotB]
        elseif inst.op == Sub
            slots[inst.out_slot] = slots[inst.in_slotA] - slots[inst.in_slotB]
        elseif inst.op == Mul
            slots[inst.out_slot] = slots[inst.in_slotA] * slots[inst.in_slotB]
        elseif inst.op == Div
            slots[inst.out_slot] = slots[inst.in_slotA] / slots[inst.in_slotB]
        elseif inst.op == Min
            slots[inst.out_slot] = min(slots[inst.in_slotA], slots[inst.in_slotB])
        elseif inst.op == Max
            slots[inst.out_slot] = max(slots[inst.in_slotA], slots[inst.in_slotB])
        else
            error("run_gpu -- unhandled op")
        end
    end
    
    if slots[1] <= 0.
        @inbounds voxels_d[x, y, z] = true
    else 
        @inbounds voxels_d[x, y, z] = false
    end

    return
end

function voxelize_gpu(dim :: Int, tape :: Tape, grid_pos :: Vec3{Float64}, grid_size :: Float64)
    @assert tape.slot_count <= MAX_SLOT_COUNT
    @assert 0 < dim <= 1024
    
    instructions_d = CuArray(tape.instructions)
    constants_d = CuArray(tape.constant_pool)
    voxels_d = CuArray{Bool}(undef, (dim, dim, dim))

    @show t = min(dim, 8)
    @show b = div(dim, t)
    @assert dim == b * t
    CUDA.@sync begin 
        @cuda blocks=(b, b, b) threads=(t, t, t) voxelize_kernel(
            dim, 
            grid_pos,
            grid_size,
            instructions_d, 
            constants_d, 
            voxels_d)
    end

    return Array(voxels_d)
end

function main()
    # Create the shape.
    shape = menger_sponge(2)
    shape = Shapes.scale(shape, 12.)
    shape = Shapes.rotateX(shape, pi / 4)
    shape = Shapes.rotateY(shape, pi / 4)
    shape = Shapes.rotateZ(shape, pi / 4)
    shape &= -Shapes.sphere(Vec3(5.0, 0.0, 6.0), 7.0)
    
    # Voxelize
    tape = node_to_tape(Nodes.constant_fold(shape))
    
    dim = 256
    voxels_raw = voxelize_gpu(dim, tape, Vec3s.full(-10.), 20.)

    voxels = Voxels.empty(Vec3s.full(-10.), 20., dim)
    for x in 1:dim 
        for y in 1:dim
            for z in 1:dim
                voxels.mask[x, y, z] = voxels_raw[x, y, z]
            end
        end
    end

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
    render!(img, camera, voxels, tape)

    # Display the image
    imshow(img)
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