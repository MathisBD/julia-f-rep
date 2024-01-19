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

function translate_aux(s :: Node, ofs :: Vec3{Node}) :: Node
    return s(axes - ofs)
end
function translate(s :: T1, ofs :: Vec3{T2}) :: Node where {T1 <: Union{Float64, Node}, T2 <: Union{Float64, Node}}
    return translate_aux(convert(Node, s), convert(Vec3{Node}, ofs))
end

function scale_aux(s :: Node, factor :: Vec3{Node}) :: Node
    return s(axes / factor)
end
function scale_aux(s :: Node, factor :: Node) :: Node
    return s(axes / factor)
end
function scale(s :: T1, factor :: Vec3{T2}) :: Node where {T1 <: Union{Float64, Node}, T2 <: Union{Float64, Node}}
    return scale_aux(convert(Node, s), convert(Vec3{Node}, factor))
end
function scale(s :: T1, factor :: T2) :: Node where {T1 <: Union{Float64, Node}, T2 <: Union{Float64, Node}}
    return scale_aux(convert(Node, s), convert(Node, factor))
end

end
