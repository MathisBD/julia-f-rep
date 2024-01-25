# A library of basic shapes and csg operators,
# represented as Nodes.
module Shapes
using ..Vec3s, ..Nodes

empty = Node(1.0)
 
axes = Vec3{Node}(Node(X), Node(Y), Node(Z))

# Shape intersection and union.
Base.:&(a :: Node, b :: Node) = max(a, b)
Base.:|(a :: Node, b :: Node) = min(a, b)


function sphere_aux(center :: Vec3{Node}, radius :: Node) :: Node
    return dist2(center, axes) - radius^2
end
function sphere(center :: Vec3{T1}, radius :: T2) :: Node where {T1 <: Union{Float64, Node}, T2 <: Union{Float64, Node}}
    return sphere_aux(convert(Vec3{Node}, center), convert(Node, radius))
end

function cube_aux(low_vertex :: Vec3{Node}, size :: Node) :: Node 
    let half_size = size / 2.0
        dx = (Node(X) - low_vertex.x - half_size)^2 - half_size^2
        dy = (Node(Y) - low_vertex.y - half_size)^2 - half_size^2
        dz = (Node(Z) - low_vertex.z - half_size)^2 - half_size^2
        return dx & dy & dz
    end
end
function cube(low_vertex :: Vec3{T1}, size :: T2) :: Node where {T1 <: Union{Float64, Node}, T2 <: Union{Float64, Node}}
    return cube_aux(convert(Vec3{Node}, low_vertex), convert(Node, size))
end

function translate(shape :: Node, ofs :: Vec3{Node}) :: Node
    return shape(axes - ofs)
end
function translate(shape :: Node, ofs :: Vec3{Float64}) :: Node
    return translate(shape, convert(Vec3{Node}, ofs))
end

function scale(shape :: Node, factor :: Vec3{Node}) :: Node
    return shape(axes / factor)
end
function scale(shape :: Node, factor :: Vec3{Float64}) :: Node
    return scale(shape, convert(Vec3{Node}, factor))
end

function scale(shape :: Node, factor :: Node) :: Node
    return shape(axes / factor)
end
function scale(shape :: Node, factor :: Float64) :: Node
    return scale(shape, convert(Node, factor))
end

function rotateX(shape :: Node, angle :: Node) :: Node
    c = cos(-angle)
    s = sin(-angle)
    return shape(Vec3(Node(X), c * Node(Y) - s * Node(Z), s * Node(Y) + c * Node(Z)))
end
function rotateX(s :: Node, angle :: Float64) :: Node
    return rotateX(s, convert(Node, angle))
end

function rotateY(shape :: Node, angle :: Node) :: Node
    c = cos(-angle)
    s = sin(-angle)
    return shape(Vec3(c * Node(X) + s * Node(Z), Node(Y), -s * Node(X) + c * Node(Z)))
end
function rotateY(s :: Node, angle :: Float64) :: Node
    return rotateY(s, convert(Node, angle))
end

function rotateZ(shape :: Node, angle :: Node) :: Node
    c = cos(-angle)
    s = sin(-angle)
    return shape(Vec3(c * Node(X) - s * Node(Y), s * Node(X) + c * Node(Y), Node(Z)))
end
function rotateZ(s :: Node, angle :: Float64) :: Node
    return rotateZ(s, convert(Node, angle))
end
    

end
