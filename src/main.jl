include(joinpath(dirname(Base.active_project()), "src", "utils", "vec.jl"))
include(joinpath(dirname(Base.active_project()), "src", "utils", "dot_graph.jl"))

include(joinpath(dirname(Base.active_project()), "src", "csg", "csg.jl"))
include(joinpath(dirname(Base.active_project()), "src", "voxelizer", "voxelizer.jl"))
include(joinpath(dirname(Base.active_project()), "src", "raytracer", "raytracer.jl"))

using ImageView, Images
using .Vec3s, .Voxels, .Cameras, .NaiveRaytracer, .NaiveVoxelizer


# Create the grid.
voxels = Voxels.empty(Vec3(-10.0, -10.0, -10.0), 20.0, 128)

# Voxelize
shape = Shapes.sphere(Vec3s.full(-10.0), 15.0)
tape = Tapes.node_to_tape(shape)
voxelize!(tape, voxels)

# Create a target image.
width = 200
height = 150
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
render!(img, camera, voxels)

# Display the image
imshow(img)

function test()
    for i in 1:1000
        render!(img, camera, voxels)
    end
end

# render! (200*150 pixels), (64^3 voxels), sphere x^2 + y^2 + z^2 <= 45^2
# --> 9.9ms


# Create a Menger sponge of depth n>=0,
# with side length 1 and centered at the origin
# function menger_sponge(n)
#    @assert n >= 0
#    if n == 0
#        return Shapes.cube(Vec3(0.0, 0.0, 0.0), 1.0)
#    else 
#        m = menger_sponge(n-1)
#        res = Shapes.empty
#        for dx in -1.0:1.0
#            for dy in -1.0:1.0
#                for dz in -1.0:1.0
#                    if abs(dx) + abs(dy) + abs(dz) > 1.0
#                        res |= Shapes.translate(m, Vec3(dx, dy, dz))
#                    end
#                end
#            end
#        end
#        return Shapes.scale(res, 1.0/3.0)
#    end
#end
#
# menger sponge 3 @ 0.1, 0.1, 0.1 ==> 315ms +- 64ms