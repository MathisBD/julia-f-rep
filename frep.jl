include("vec.jl")
include("dot_graph.jl")

module Nodes
using ..Vec3s, ..DotGraph

export Vec3s
export Op, Node, X, Y, Z, Const, Add, Sub, Mul, Min, Max
export topoiter, topomap, visualize

@enum Op X Y Z Const Add Sub Mul Min Max

struct Node
   op :: Op
   inputs :: Vector{Node}
   constant :: Number

   Node(c :: Number) = new(Const, [], c)
   Node(op :: Op) = new(op, [], 0)
   Node(op :: Op, inp1 :: Node) = new(op, [inp1], 0)
   Node(op :: Op, inp1 :: Node, inp2 :: Node) = new(op, [inp1, inp2], 0)
end

function topoiter(f :: Function, root :: Node) 
    visited = Set() 
    
    function dfs(node)
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

# T is the type of the result of f.
function topomap(T :: Type, f :: Function, root :: Node) :: T
    result = Dict{Node, T}()
    
    function helper(node :: Node)
        input_results = [result[child] for child in node.inputs]
        result[node] = f(node, input_results)
    end
    
    topoiter(helper, root)
    return result[root]
end

# Render the dot graph corresponding to a Node.
function visualize(root :: Node)
    graph = DotGraph.Graph(true)
    
    function helper(node :: Node, children_ids :: Vector{Int})
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

Base.:+(a :: Node, b :: Node) = Node(Add, a, b)
Base.:-(a :: Node, b :: Node) = Node(Sub, a, b)
Base.:*(a :: Node, b :: Node) = Node(Mul, a, b)
Base.min(a :: Node, b :: Node) = Node(Min, a, b)
Base.max(a :: Node, b :: Node) = Node(Max, a, b)
function Base.:^(a :: Node, n :: Int)
    @assert n >= 0
    if n == 0
        return Node(1)
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

end

#Base.:-(a :: Node, b :: Node) = Sub(Ref{Node}(a), Ref{Node}(b))
#Base.:*(a :: Node, b :: Node) = Mul(Ref{Node}(a), Ref{Node}(b))
#Base.min(a :: Node, b :: Node) = Min(Ref{Node}(a), Ref{Node}(b))
#Base.max(a :: Node, b :: Node) = Max(Ref{Node}(a), Ref{Node}(b))
#function Base.:^(a :: Node, n :: Int)
#    @assert n >= 0
#    if n == 0
#        return Const(1)
#    elseif n == 1
#        return a
#    elseif n % 2 == 0
#        b = a^(div(n, 2))
#        return Mul(Ref{Node}(b), Ref{Node}(b))
#    else
#        b = a^(div(n, 2))
#        return a * Mul(Ref{Node}(b), Ref{Node}(b))
#    end
#end
#
## Formally replace the arguments X(), Y() and Z() of a node.
#function (n::Node)(v::Vec3{Node})
#    @match n begin
#        X() => v.x
#        Y() => v.y
#        Z() => v.z
#        Const(c) => Const(c)
#        Add(n1, n2) => Add(n1[]())
#        Sub(n1, n2) => 
#        Mul(n1, n2) => 
#        Min(n1, n2) => 
#        Max(n1, n2) => 
#    end
#end


## Evaluate a node at a point.
#function eval_node(n::Node, v::Vec3)
#    @match n begin
#        X() => v.x
#        Y() => v.y
#        Z() => v.z
#        Const(c) => c
#        Add(n1, n2) => n1[](v) + n2[](v)
#        Sub(n1, n2) => n1[](v) - n2[](v)
#        Mul(n1, n2) => n1[](v) * n2[](v)
#        Min(n1, n2) => min(n1[](v), n2[](v))
#        Max(n1, n2) => max(n1[](v), n2[](v))
#    end
#end

# A library of basic shapes and csg operators,
# represented as Nodes.
#module Shapes
#using ..Vec3s, ..Nodes
#export Vec3s
#
#axes = Vec3{Node}(X(), Y(), Z())
#
## Shape intersection.
#Base.:&(a :: Node, b :: Node) = max(a, b)
## Shape union.
#Base.:|(a :: Node, b :: Node) = min(a, b)
#
#sphere(center :: Vec3{Node}, radius :: Node) = dist2(center, axes) - radius^2
#cube(low_vertex :: Vec3{Node}, size :: Node) = 
#    let dx = (X() - low_vertex.x - size / Const(2))^2 - (size / Const(2))^2
#        dy = (Y() - low_vertex.y - size / Const(2))^2 - (size / Const(2))^2
#        dz = (Z() - low_vertex.z - size / Const(2))^2 - (size / Const(2))^2
#        dx & dy & dz
#    end
#
#translate(s :: Node, ofs :: Vec3{Node}) = s()
#
#end
#
#
## Create a Menger sponge of depth n>=0,
## with side length 1 and centered at the origin
#function menger_sponge(n)
#    @assert n >= 0
#    if n == 0
#        Shapes.cube(Vec3(0.0, 0.0, 0.0), 1.0)
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
#
### Represent functions from T^3 to T as tapes.
##module FRepTape
### TODO : define the [Tape] struct
##end
##
### TEST : create the Menger sponge. 
### Benchmark how long it takes to evaluate at a few different points,
### both for the basic version and the Tree version (naive and simplified), 
### and finally the Tape version. 