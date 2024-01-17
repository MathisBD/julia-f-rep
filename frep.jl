include("vec.jl")
include("dot_graph.jl")

module Nodes
using ..Vec3s, ..DotGraph
import Base: convert

export Op, Node, X, Y, Z, Const, Add, Sub, Mul, Div, Min, Max
export topoiter, topomap, visualize

@enum Op X Y Z Const Add Sub Mul Div Min Max

struct Node{T <: Number}
    op :: Op
    inputs :: Vector{Node{T}}
    constant :: T
 
    Node{T}(c :: T) where {T <: Number} = 
        new{T}(Const, [], c)
    Node{T}(op :: Op) where {T <: Number} = 
        new{T}(op, [], zero(T))
    Node{T}(op :: Op, inp1 :: Node{T}) where {T <: Number} = 
        new{T}(op, [inp1], zero(T))
    Node{T}(op :: Op, inp1 :: Node{T}, inp2 :: Node{T}) where {T <: Number} = 
        new{T}(op, [inp1, inp2], zero(T))
end

# Convert a number to a node.
Base.convert(:: Type{Node{T}}, x :: T) where {T <: Number} = Node{T}(x)

Base.:+(a :: Node{T}, b :: Node{T}) where {T <: Number} = Node{T}(Add, a, b)
Base.:+(a :: Node{T}, b :: T) where {T <: Number} = Node{T}(Add, a, convert(Node{T}, b))
Base.:+(a :: T, b :: Node{T}) where {T <: Number} = Node{T}(Add, convert(Node{T}, a), b)

Base.:-(a :: Node{T}, b :: Node{T}) where {T <: Number} = Node{T}(Sub, a, b)
Base.:-(a :: Node{T}, b :: T) where {T <: Number} = Node{T}(Sub, a, convert(Node, b))
Base.:-(a :: T, b :: Node{T}) where {T <: Number} = Node{T}(Sub, convert(Node{T}, a), b)

Base.:*(a :: Node{T}, b :: Node{T}) where {T <: Number} = Node{T}(Mul, a, b)
Base.:*(a :: Node{T}, b :: T) where {T <: Number} = Node{T}(Mul, a, convert(Node{T}, b))
Base.:*(a :: T, b :: Node{T}) where {T <: Number} = Node{T}(Mul, convert(Node{T}, a), b)

Base.:/(a :: Node{T}, b :: Node{T}) where {T <: Number} = Node{T}(Div, a, b)
Base.:/(a :: Node{T}, b :: T) where {T <: Number} = Node{T}(Div, a, convert(Node{T}, b))
Base.:/(a :: T, b :: Node{T}) where {T <: Number} = Node{T}(Div, convert(Node{T}, a), b)

Base.min(a :: Node{T}, b :: Node{T}) where {T <: Number} = Node{T}(Min, a, b)
Base.min(a :: Node{T}, b :: T) where {T <: Number} = Node{T}(Min, a, convert(Node{T}, b))
Base.min(a :: T, b :: Node{T}) where {T <: Number} = Node{T}(Min, convert(Node{T}, a), b)

Base.max(a :: Node{T}, b :: Node{T}) where {T <: Number} = Node{T}(Max, a, b)
Base.max(a :: Node{T}, b :: T) where {T <: Number} = Node{T}(Max, a, convert(Node{T}, b))
Base.max(a :: T, b :: Node{T}) where {T <: Number} = Node{T}(Max, convert(Node{T}, a), b)

function Base.:^(a :: Node{T}, n :: Int) where {T <: Number} 
    @assert n >= 0
    if n == 0
        return Node{T}(one(T))
    elseif n == 1
        return a
    elseif n % 2 == 0
        b = a^(div(n, 2))
        return b * b
    else
        b = a^(div(n, 2))
        return a * (b * b)
    end
end

# Apply a function f on root and all of its subnodes.
function topoiter(f :: Function, root :: Node{T}) where {T} 
    visited = Set{Node{T}}() 
    
    function dfs(node :: Node{T})
        if !(node in visited)
            push!(visited, node)
            for child in node.inputs
                dfs(child)
            end
            f(node)
        end
    end

    dfs(root)
end

# Tres is the type of the result of f.
function topomap(:: Type{Tres}, f :: Function, root :: Node{Tnode}) :: Tres where {Tnode, Tres}
    result = Dict{Node{Tnode}, Tres}()
    
    function helper(node :: Node{Tnode})
        input_results = [result[child] for child in node.inputs]
        result[node] = f(node, input_results)
    end
    
    topoiter(helper, root)
    return result[root]
end

# Render the dot graph corresponding to a Node.
function visualize(root :: Node{T}) where {T}
    graph = DotGraph.Graph(true)
    
    function helper(node :: Node{T}, children_ids :: Vector{Int})
        # Compute the label of the node.
        label = node == root ? "ROOT-" : ""
        if node.op == Const
            label *= "$(node.constant)"
        else
            label *= "$(node.op)"
        end

        # Add a node to the graph.
        id = DotGraph.add_node!(graph, label)

        # Add the edges to the child nodes.
        for child_id in children_ids 
            DotGraph.add_edge!(graph, id, child_id)
        end

        return id
    end

    topomap(Int, helper, root)
    return DotGraph.render(graph)
end

# Formally replace the arguments X, Y and Z of a node.
# This can be used to replace X, Y and Z by other nodes,
# or by numeric values to compute the numeric value of the node.
function (root :: Node{Tnode})(v :: Vec3{T}) :: T where {Tnode, T}
    function helper(n :: Node{Tnode}, children :: Vector{T}) :: T
        if n.op == X
            return v.x
        elseif n.op == Y
            return v.y
        elseif n.op == Z
            return v.z
        elseif n.op == Const
            return convert(T, n.constant)
        elseif n.op == Add
            return children[1] + children[2]
        elseif n.op == Sub
            return children[1] - children[2]
        elseif n.op == Mul
            return children[1] * children[2]
        elseif n.op == Div
            return children[1] / children[2]
        elseif n.op == Min
            return min(children[1], children[2])
        elseif n.op == Max
            return max(children[1], children[2])
        else
            error("unhandled value $(n.op) :: $(typeof(n.op))")
        end
    end

    return topomap(T, helper, root)
end

end


# A library of basic shapes and csg operators,
# represented as Nodes.
module Shapes
using ..Vec3s, ..Nodes

empty(::Type{T}) where {T <: Number} = Node{T}(one(T))

axes(::Type{T}) where {T <: Number} = Vec3{Node{T}}(Node{T}(X), Node{T}(Y), Node{T}(Z))

# Shape intersection and union.
Base.:&(a :: Node{T}, b :: Node{T}) where {T} = max(a, b)
Base.:|(a :: Node{T}, b :: Node{T}) where {T} = min(a, b)


function sphere(center :: Vec3{Node{T}}, radius :: Node{T}) :: Node{T} where {T}
    return dist2(center, axes(T)) - radius^2
end
#function sphere(center :: Vec3{T}, radius :: Node{T}) :: Node{T} where {T <: Number}
#    return sphere(convert(Vec3{Node{T}}, center), radius)
#end

function cube_aux(low_vertex :: Vec3{Node{T}}, size :: Node{T}) :: Node{T} where {T} 
    let half_size = size / convert(T, 2)
        dx = (Node{T}(X) - low_vertex.x - half_size)^2 - half_size^2
        dy = (Node{T}(Y) - low_vertex.y - half_size)^2 - half_size^2
        dz = (Node{T}(Z) - low_vertex.z - half_size)^2 - half_size^2
        return dx & dy & dz
    end
end
function cube(low_vertex :: Vec3{T1}, size :: T2) :: Node{T} where {T <: Number, T1 <: Union{T, Node{T}}, T2 <: Union{T, Node{T}}}
    return cube_aux(convert(Vec3{Node{T}}, low_vertex), convert(Node{T}, size))
end

function translate_aux(s :: Node{T}, ofs :: Vec3{Node{T}}) :: Node{T} where {T}
    return s(axes(T) - ofs)
end
function translate(s :: T1, ofs :: Vec3{T2}) :: Node{T} where {T <: Number, T1 <: Union{T, Node{T}}, T2 <: Union{T, Node{T}}}
    return translate_aux(convert(Node{T}, s), convert(Vec3{Node{T}}, ofs))
end

function scale_aux(s :: Node{T}, factor :: Vec3{Node{T}}) :: Node{T} where {T}
    return s(axes(T) / factor)
end
function scale_aux(s :: Node{T}, factor :: Node{T}) :: Node{T} where {T}
    return s(axes(T) / factor)
end
#function scale(s :: T1, factor :: Vec3{T2}) :: Node where {T1 <: Union{Node, Number}, T2 <: Union{Node, Number}}
#    return scale(convert(Node, s), convert(Vec3{Node}, factor))
#end
#function scale(s :: T1, factor :: T2) :: Node where {T1 <: Union{Node, Number}, T2 <: Union{Node, Number}}
#    return scale(convert(Node, s), Vec3(convert(Node, factor), convert(Node, factor), convert(Node, factor)))
#end

end

using .Shapes, .Vec3s, .Nodes

# Create a Menger sponge of depth n>=0,
# with side length 1 and centered at the origin
function menger_sponge(n)
    @assert n >= 0
    if n == 0
        return Shapes.cube(Vec3(0.0, 0.0, 0.0), Node{Float64}(1.0))
    else 
        m = menger_sponge(n-1)
        res = Shapes.empty(Float64)
        for dx in -1.0:1.0
            for dy in -1.0:1.0
                for dz in -1.0:1.0
                    if abs(dx) + abs(dy) + abs(dz) > 1.0
                        res |= Shapes.translate(m, Vec3(dx, dy, dz))
                    end
                end
            end
        end
        return Shapes.scale(res, Node{Float64}(1.0/3.0))
    end
end
