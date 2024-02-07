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

const MAX_SLOT_COUNT = 32


function run_tape(
    id,
    instructions_d :: CuDeviceArray{Instruction}, 
    constants_d :: CuDeviceArray{Float32},
    wx :: Float32, 
    wy :: Float32, 
    wz :: Float32) :: Float32

    slots = zeros(MVector{MAX_SLOT_COUNT, Float32})
    slots[1] = wx
    slots[2] = wy    
    slots[3] = wz
    
    # Execute the tape instructions.
    @inbounds for i in 1:length(instructions_d)
        inst = instructions_d[i]
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
    return slots[1]
end

function voxelize_kernel(
    dim :: Int,
    grid_pos :: Vec3{Float32},
    grid_size :: Float32,
    instructions_d :: CuDeviceArray{Instruction}, 
    constants_d :: CuDeviceArray{Float32},
    voxels_d :: CuDeviceArray{Bool, 3})

    # id in range [0, dim^3 - 1]
    id = (blockIdx().x - 1) * blockDim().x + threadIdx().x - 1

    # x, y, z in range [0, dim - 1]
    x = id % dim
    y = div(id, dim) % dim
    z = div(id, dim^2)

    wx = x * grid_size / dim + grid_pos.x 
    wy = y * grid_size / dim + grid_pos.y
    wz = z * grid_size / dim + grid_pos.z
    res = run_tape(id, instructions_d, constants_d, wx, wy, wz)

    if res <= 0f0
        @inbounds voxels_d[x+1, y+1, z+1] = true
    else 
        @inbounds voxels_d[x+1, y+1, z+1] = false
    end

    return
end

function voxelize_gpu(dim :: Int, tape :: Tape, grid_pos :: Vec3{Float32}, grid_size :: Float32)
    @assert tape.slot_count <= MAX_SLOT_COUNT
    @assert 0 < dim <= 2048
    
    instructions_d = CuArray(tape.instructions)
    constants_d = CuArray(Float32.(tape.constant_pool))
    voxels_d = CuArray{Bool}(undef, (dim, dim, dim))

    # Compile the kernel.
    kernel = @cuda launch=false voxelize_kernel(
        dim, grid_pos, grid_size, instructions_d, constants_d, voxels_d)

    # Compute the launch configuration.
    @show threads = min(dim^3, 256)
    @show blocks = div(dim^3, threads)
    @assert dim^3 == blocks * threads

    # Run the kernel.
    CUDA.@sync begin 
        kernel(dim, grid_pos, grid_size, instructions_d, constants_d, voxels_d; 
            threads, blocks)
    end

    return Array(voxels_d)
end

function benchmark_voxelize_gpu(dim :: Int, tape :: Tape, grid_pos :: Vec3{Float32}, grid_size :: Float32)
    @assert tape.slot_count <= MAX_SLOT_COUNT
    @assert 0 < dim <= 2048
    
    instructions_d = CuArray(tape.instructions)
    constants_d = CuArray(Float32.(tape.constant_pool))
    voxels_d = CuArray{Bool}(undef, (dim, dim, dim))

    # Compile the kernel.
    kernel = @cuda launch=false voxelize_kernel(
        dim, grid_pos, grid_size, instructions_d, constants_d, voxels_d)

    # Compute the launch configuration.
    @show threads = min(dim^3, 256)
    @show blocks = div(dim^3, threads)
    @assert dim^3 == blocks * threads

    function bench()
        # Run the kernel.
        CUDA.@sync begin 
            kernel(dim, grid_pos, grid_size, instructions_d, constants_d, voxels_d; 
                threads, blocks)
        end
    end

    @benchmark $bench()
end


#function main()
    # Create the shape.
    shape = menger_sponge(1)
    shape = Shapes.scale(shape, 12.)
    shape = Shapes.rotateX(shape, pi / 4)
    shape = Shapes.rotateY(shape, pi / 4)
    shape = Shapes.rotateZ(shape, pi / 4)
    shape &= -Shapes.sphere(Vec3(5.0, 0.0, 6.0), 7.0)
    
    # Voxelize
    tape = node_to_tape(Nodes.constant_fold(shape))
    
    dim = 256
    grid_pos = Vec3s.full(-10f0)
    grid_size = 20f0
    voxels_raw = voxelize_gpu(dim, tape, grid_pos, grid_size)

    # Copy the voxels
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
#end
