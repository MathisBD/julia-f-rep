# A representation of expression trees that is more efficient to evaluate.
module Tapes
using ..Vec3s, Printf
import ..Nodes
import Base: print

export Op, Copy, LoadConst, Sin, Cos, Exp, Neg, Sqrt, Add, Sub, Mul, Div, Min, Max, SMin, SMax
export Instruction, Tape, node_to_tape, run


@enum Op::UInt8 Copy LoadConst Sin Cos Exp Neg Sqrt Add Sub Mul Div Min Max SMin SMax

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
    constant_pool :: Vector{Float64}
    slot_count :: Int
end

function node_op_to_tape_op(op :: Nodes.Op)
    if op == Nodes.Const    return LoadConst
    elseif op == Nodes.Neg  return Neg
    elseif op == Nodes.Sin  return Sin
    elseif op == Nodes.Cos  return Cos
    elseif op == Nodes.Exp  return Exp
    elseif op == Nodes.Sqrt return Sqrt   
    elseif op == Nodes.Add  return Add
    elseif op == Nodes.Sub  return Sub
    elseif op == Nodes.Mul  return Mul
    elseif op == Nodes.Div  return Div
    elseif op == Nodes.Min  return Min
    elseif op == Nodes.Max  return Max
    elseif op == Nodes.SMin  return SMin
    elseif op == Nodes.SMax  return SMax
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
    last_slot = -1
    for (i, n) in enumerate(nodes)
        local op :: Op
        local in_slotA :: UInt8 = 0
        local in_slotB :: UInt8 = 0
        local out_slot :: UInt8 = 0

        # Axis ops : these don't generate any instruction.             
        if n.op == Nodes.X 
            last_slot = 1
            continue
        elseif n.op == Nodes.Y
            last_slot = 2
            continue
        elseif n.op == Nodes.Z
            last_slot = 3
            continue
        # LoadConst
        elseif n.op == Nodes.Const
            @assert 1 <= constant_idx[n.constant] <= length(constant_pool)
            in_slotA = constant_idx[n.constant]
        # Nodes with one input.
        elseif Nodes.has_one_input(n.op)
            in_slotA = get_current_slot(n.inputs[1])
            
            @assert liveness[node_idx[n.inputs[1]]] >= i
            if liveness[node_idx[n.inputs[1]]] == i
                slots[in_slotA] = nothing
            end            
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
        last_slot = out_slot
    end

    # In some cases (for instance root == Node(Y)), there is no instruction
    # generated and the output is not in slot 1.
    @assert last_slot >= 1
    if last_slot != 1
        push!(instructions, Instruction(Copy, 1, last_slot, 0))
    end

    @assert 1 <= length(slots)
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

# The type T we use to run the computations must be a number type
# that supports conversion from Float64.
function run(tape :: Tape, x :: T, y :: T, z :: T) :: T where {T <: Number}
    # Initially fill the slots with a dummy value.
    slots = fill(x, tape.slot_count)
    slots[1] = x
    slots[2] = y
    slots[3] = z
    
    for inst in tape.instructions
        if inst.op == Copy
            slots[inst.out_slot] = slots[inst.in_slotA]
        elseif inst.op == LoadConst
            slots[inst.out_slot] = convert(T, tape.constant_pool[inst.in_slotA])
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
        elseif inst.op == SMin
            slots[inst.out_slot] = smooth_min(slots[inst.in_slotA], slots[inst.in_slotB])
        elseif inst.op == SMax
            slots[inst.out_slot] = smooth_max(slots[inst.in_slotA], slots[inst.in_slotB])
        else
            error("Unhandled op : $(inst.op) :: $(typeof(inst.op))")
        end
    end

    return slots[1]
end

end
