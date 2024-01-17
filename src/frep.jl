include("vec.jl")
include("dot_graph.jl")

module Nodes
using ..Vec3s, ..DotGraph
import Base: convert, show

export Op, Node, X, Y, Z, Const, Add, Sub, Mul, Div, Min, Max
export topoiter, topomap, visualize, constant_fold

@enum Op X Y Z Const Add Sub Mul Div Min Max

struct Node
    op :: Op
    inputs :: Vector{Node}
    constant :: Float64
 
    Node(c :: Float64) = new(Const, [], c)
    Node(op :: Op) = new(op, [], 0.0)
    Node(op :: Op, inp1 :: Node) = new(op, [inp1], 0.0)
    Node(op :: Op, inp1 :: Node, inp2 :: Node) = new(op, [inp1, inp2], 0.0)
end

function Base.show(io :: IO, n :: Node)
    if n.op == X || n.op == Y || n.op == Z
        print(io, "Node($(n.op))")
    elseif n.op == Const
        print(io, "Node($(n.constant))")
    else 
        print(io, "Node($(n.op), $(n.inputs[1]), $(n.inputs[2]))")
    end
end

function has_inputs(op :: Op) :: Bool
    return op == Add || op == Sub || op == Mul || op == Div || op == Min || op == Max
end

function is_constant(node :: Node, constant :: Float64)
    return node.op == Const && node.constant == constant
end

# Apply an operator to numeric arguments. 
# This doesn't work for all operators.
function apply_op(::Val{Add}, a :: T, b :: T) :: T where {T} 
    return a + b
end
function apply_op(::Val{Sub}, a :: T, b :: T) :: T where {T} 
    return a - b
end
function apply_op(::Val{Mul}, a :: T, b :: T) :: T where {T} 
    return a * b
end
function apply_op(::Val{Div}, a :: T, b :: T) :: T where {T} 
    return a / b
end
function apply_op(::Val{Min}, a :: T, b :: T) :: T where {T} 
    return min(a, b)
end
function apply_op(::Val{Max}, a :: T, b :: T) :: T where {T} 
    return max(a, b)
end

# Convert a float to a node.
convert(::Type{Node}, x :: Float64) = Node(x)

Base.:+(a :: Node, b :: Node) = Node(Add, a, b)
Base.:+(a :: Node, b :: Float64) = Node(Add, a, Node(b))
Base.:+(a :: Float64, b :: Node) = Node(Add, Node(a), b)

Base.:-(a :: Node, b :: Node) = Node(Sub, a, b)
Base.:-(a :: Node, b :: Float64) = Node(Sub, a, Node(b))
Base.:-(a :: Float64, b :: Node) = Node(Sub, Node(a), b)

Base.:*(a :: Node, b :: Node) = Node(Mul, a, b)
Base.:*(a :: Node, b :: Float64) = Node(Mul, a, Node(b))
Base.:*(a :: Float64, b :: Node) = Node(Mul, Node(a), b)

Base.:/(a :: Node, b :: Node) = Node(Div, a, b)
Base.:/(a :: Node, b :: Float64) = Node(Div, a, Node(b))
Base.:/(a :: Float64, b :: Node) = Node(Div, Node(a), b)

Base.min(a :: Node, b :: Node) = Node(Min, a, b)
Base.min(a :: Node, b :: Float64) = Node(Min, a, Node(b))
Base.min(a :: Float64, b :: Node) = Node(Min, Node(a), b)

Base.max(a :: Node, b :: Node) = Node(Max, a, b)
Base.max(a :: Node, b :: Float64) = Node(Max, a, Node(b))
Base.max(a :: Float64, b :: Node) = Node(Max, Node(a), b)

function Base.:^(a :: Node, n :: Int) 
    @assert n >= 0
    if n == 0
        return Node(1.0)
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
function topoiter(f :: Function, root :: Node)
    visited = Set{Node}() 
    
    function dfs(node :: Node)
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
function topomap(:: Type{T}, f :: Function, root :: Node) :: T where {T}
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

# Formally replace the arguments X, Y and Z of a node.
# This can be used to replace X, Y and Z by other nodes,
# or by numeric values to compute the numeric value of the node.
# T is either Float64 or Node.
function (root :: Node)(v :: Vec3{T}) :: T where {T <: Union{Float64, Node}}
    function helper(n :: Node, children :: Vector{T}) :: T
        if n.op == X
            return v.x
        elseif n.op == Y
            return v.y
        elseif n.op == Z
            return v.z
        elseif n.op == Const
            return convert(T, n.constant)
        else
            return apply_op(Val(n.op), children...)
        end
    end

    return topomap(T, helper, root)
end

include("frep_simplify.jl")

end


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


# A representation of expression trees that is more efficient to evaluate.
module Tapes
using ..Nodes, ..Vec3s

@enum Opcode::UInt8 Copy LoadConst Sin Cos Exp Neg Sqrt Add Sub Mul Div Min Max

struct Instruction 
    opcode :: Opcode
    out_slot :: UInt8
    in_slotA :: UInt8
    in_slotB :: UInt8
end

struct Tape
    instructions :: Vector{Instruction}
    constant_pool :: Vector{Float64}
    slot_count :: Int
end

function node_to_tape(root :: Node) :: Tape
    # Get a topological sort of the nodes in the tree.
    nodes = Node[]
    topoiter(n -> push!(nodes, n), root)
    # And the index of each node in this topo sort.
    node_idx = Dict{Node, Int}
    for (i, n) in enumerate(nodes)
        node_idx[node] = i
    end

    # Build the constant pool.
    constant_pool = Float64[]
    constant_idx = Dict{Node}
    for n in nodes
        if n.op == Const && !(n.constant in constant_pool)
            push!(constant_pool, n.constant)
            constant_idx[n.constant] = length(constant_pool)
        end
    end

    # Compute the liveness of every node, i.e. 
    # the index of the last node that uses its result
    # (or -1 if its result is never used).
    liveness = repeat([-1], length(nodes))
    for (i, n) in enumerate(nodes)
        if has_inputs(n)
            for inp in n.inputs 
                # The max is not strictly necessary but clarifies what we are doing.
                liveness[node_idx[inp]] = max(i, liveness[node_idx[inp]])
            end
        end
    end

    # Build the instruction for each node.
    slots = Node[] # TODO Node option

    function get_free_slot()
        for i in eachindex(slots)
            if slots[i] = 
    end

end

end

using .Shapes, .Vec3s, .Nodes

# Create a Menger sponge of depth n>=0,
# with side length 1 and centered at the origin
function menger_sponge(n)
    @assert n >= 0
    if n == 0
        return Shapes.cube(Vec3(0.0, 0.0, 0.0), 1.0)
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
        return Shapes.scale(res, 1.0/3.0)
    end
end

# menger sponge 3 @ 0.1, 0.1, 0.1 ==> 315ms +- 64ms