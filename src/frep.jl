include("vec.jl")
include("dot_graph.jl")

module Nodes
using ..Vec3s, ..DotGraph
import Base: convert, show

export Op, Node, X, Y, Z, Const, Add, Sub, Mul, Div, Min, Max
export topoiter, topomap, visualize, constant_fold, merge_axes

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

function has_one_input(op :: Op) :: Bool
    return false
end

function has_two_inputs(op :: Op) :: Bool
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
using ..Vec3s, Printf
import ..Nodes
import Base: print

@enum Op::UInt8 Copy LoadConst Sin Cos Exp Neg Sqrt Add Sub Mul Div Min Max

# An instruction should be as small as possible (32 bits at the moment).
struct Instruction 
    op :: Op
    out_slot :: UInt8
    # For a LoadConst instruction, this is the index in the constant pool
    in_slotA :: UInt8 
    in_slotB :: UInt8
end

struct Tape
    instructions :: Vector{Instruction}
    # LoadConst instructions don't contain a float value,
    # but instead contain an index in the constant pool.
    # The axis instructions (X, Y, Z) get converted to LoadConst
    # instructions : the values (x, y, z) are the first three constants in the pool.
    constant_pool :: Vector{Float64}
    slot_count :: Int
end

function node_op_to_tape_op(op :: Nodes.Op)
    if op == Nodes.X || op == Nodes.Y || op == Nodes.Z || op == Nodes.Const
        return LoadConst
    elseif op == Nodes.Add
        return Add
    elseif op == Nodes.Sub
        return Sub
    elseif op == Nodes.Mul
        return Mul
    elseif op == Nodes.Div
        return Div
    elseif op == Nodes.Min
        return Min
    elseif op == Nodes.Max
        return Max
    else
        error("unhandled $op :: $(typeof(op))")
    end
end

function node_to_tape(root :: Nodes.Node) :: Tape
    # Get a topological sort of the nodes in the tree.
    # We also need to make sure each axis (X, Y and Z) appears only once in root. 
    nodes = Nodes.Node[]
    Nodes.topoiter(n -> push!(nodes, n), Nodes.merge_axes(root))
    # And the index of each node in this topo sort.
    node_idx = Dict{Nodes.Node, Int}()
    for (i, n) in enumerate(nodes)
        node_idx[n] = i
    end

    # Build the constant pool.
    # The first values are placeholders : they will be 
    # replaced by the values of x, y, z when we evaluate the tape.
    constant_pool = Float64[]
    constant_idx = Dict{Float64, Int}()
    for n in nodes
        if n.op == Nodes.Const && !(n.constant in constant_pool)
            push!(constant_pool, n.constant)
            constant_idx[n.constant] = length(constant_pool)
        end
    end
    # Check we can index in the constant pool using UInt8.
    @assert length(constant_pool) < 256
    
    # Compute the liveness of every node, i.e. 
    # the index of the last node that uses its result
    # (or -1 if its result is never used).
    liveness = fill(-1, length(nodes))
    for (i, n) in enumerate(nodes)
        if Nodes.has_inputs(n.op)
            for inp in n.inputs 
                # The max is not strictly necessary but clarifies what we are doing.
                liveness[node_idx[inp]] = max(i, liveness[node_idx[inp]])
            end
        end
    end

    # Build the initial contents of the slots.
    slots = Union{Nothing, Nodes.Node}[]
    local node_x :: Union{Nothing, Nodes.Node} = nothing
    local node_y :: Union{Nothing, Nodes.Node} = nothing
    local node_z :: Union{Nothing, Nodes.Node} = nothing
    for n in nodes
        if n.op == Nodes.X
            @assert isnothing(node_x)
            node_x = n
        elseif n.op == Nodes.Y
            @assert isnothing(node_y)
            node_y = n
        elseif n.op == Nodes.Z
            @assert isnothing(node_z)
            node_z = n
        end
    end
    push!(slots, node_x, node_y, node_z)
    

    function get_free_slot()
        for (i, s) in enumerate(slots)
            if isnothing(s)
                return i
            end
        end
        # No free slots : add one.
        push!(slots, nothing)
        return length(slots)
    end

    function get_current_slot(n :: Nodes.Node)
        for (i, s) in enumerate(slots)
            if s == n
                return i
            end 
        end
        error("Could not find slot for node $n")
    end

    # Build the instruction for each node.
    instructions = Instruction[]
    for (i, n) in enumerate(nodes)
        local op :: Op
        local in_slotA :: UInt8 = 0
        local in_slotB :: UInt8 = 0
        local out_slot :: UInt8 = 0

        # Axis ops : these don't generate any instruction.             
        if n.op == Nodes.X || n.op == Nodes.Y || n.op == Nodes.Z
            continue
        # LoadConst
        elseif n.op == Nodes.Const
            @assert 1 <= constant_idx[n.constant] <= length(constant_pool)
            in_slotA = constant_idx[n.constant]
        # Nodes with one input.
        elseif Nodes.has_one_input(n.op)
            error("not implemented yet")
        # Nodes with two inputs.
        elseif Nodes.has_two_inputs(n.op)
            in_slotA = get_current_slot(n.inputs[1])
            in_slotB = get_current_slot(n.inputs[2])

            @assert liveness[node_idx[n.inputs[1]]] >= i
            if liveness[node_idx[n.inputs[1]]] == i
                slots[in_slotA] = nothing
            end            
            @assert liveness[node_idx[n.inputs[2]]] >= i
            if liveness[node_idx[n.inputs[2]]] == i
                slots[in_slotB] = nothing
            end
        end

        # We get a free slot for the ouput AFTER having 
        # released the input slots. This way an instruction 
        # can store its output in a slot previously used by an input.
        op = node_op_to_tape_op(n.op)
        out_slot = get_free_slot()
        slots[out_slot] = n
        push!(instructions, Instruction(op, out_slot, in_slotA, in_slotB))
    end

    # Check the output slot is 1.
    @assert 1 <= length(slots)
    @assert instructions[length(instructions)].out_slot == 1
    # Check we can index slots using UInt8.
    @assert length(slots) < 256

    return Tape(instructions, constant_pool, length(slots))
end

function Base.print(io :: IO, tape :: Tape)
    print(io, "[+] Tape instr_count=$(length(tape.instructions)) slot_count=$(tape.slot_count)\n")
    
    for (i, inst) in enumerate(tape.instructions)
        @printf(io, "    [%i] %-10s inA=%-3i inB=%-3i out=%-3i\n", 
            i, inst.op, inst.in_slotA, inst.in_slotB, inst.out_slot)
    end

    print(io, "[+] Constant pool size=$(length(tape.constant_pool))\n")
    for (i, c) in enumerate(tape.constant_pool)
        @printf(io, "    [%i] %4.2f\n", i, c)
    end
end

function run(tape :: Tape, x :: Float64, y :: Float64, z :: Float64) :: Float64
    slots = fill(0.0, tape.slot_count)
    slots[1] = x
    slots[2] = y
    slots[3] = z

    for inst in tape.instructions
        if inst.op == Copy
            slots[inst.out_slot] = slots[inst.in_slotA]
        elseif inst.op == LoadConst
            slots[inst.out_slot] = tape.constant_pool[inst.in_slotA]
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
            error("Unhandled op : $(inst.op) :: $(typeof(inst.op))")
        end
    end

    return slots[1]
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