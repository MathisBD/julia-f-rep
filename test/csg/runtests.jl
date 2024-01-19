
include(joinpath(dirname(Base.active_project()), "src", "csg", "csg.jl"))

using .Vec3s, .Nodes, .Shapes
import .Tapes

using JCheck
import JCheck: generate, shrinkable, shrink

# FRep generators

# Return a random float64 in a reasonable range (contrary to what generate does).
function rand_float64(rng)
    return (rand(rng, Float64) - 0.5) * 10.0
end
function rand_float64(rng, n :: Int)
    return (rand(rng, Float64, n) .- 0.5) .* 10.0
end

# size is the total number of non-leaf nodes (in the tree) we should generate.
function rand_node_sized(rng, size :: Int)
    @assert size >= 0

    # Choose an operator
    local op :: Op
    if size == 0
        op = rand(rng, (X, Y, Z, Const))
    else 
        # TODO : maybe add support for partial operations (Div, Sqrt, ...)
        op = rand(rng, (Neg, Sin, Cos, Add, Sub, Mul, Min, Max))
    end

    # Axis node.
    if op == X || op == Y || op == Z 
        @assert size == 0
        return Node(op)
    # Constant node.
    elseif op == Const
        @assert size == 0
        return Node(rand_float64(rng))
    # Node with one input.
    elseif Nodes.has_one_input(op)
        @assert size > 0
        child = rand_node_sized(rng, size - 1)
        return Node(op, child)
    # Node with two inputs.
    elseif Nodes.has_two_inputs(op)
        @assert size > 0
        size_left = rand(rng, 0:(size-1))
        size_right = size - 1 - size_left
        
        left = rand_node_sized(rng, size_left)
        right = rand_node_sized(rng, size_right)
        return Node(op, left, right) 
    else 
        error("Unhandled op $(op) :: $(typeof(op))")
    end
end

function generate(rng, ::Type{Node}, n :: Int)
    max_size = 100
    return [rand_node_sized(rng, s) for s in rand(rng, 0:max_size, n)]
end

function generate(rng, ::Type{Vec3{Float64}}, n :: Int)
    xs = rand_float64(rng, n)
    ys = rand_float64(rng, n)
    zs = rand_float64(rng, n)
    return [Vec3(xyz...) for xyz in zip(xs, ys, zs)]
end

# FRep shrinking

function shrinkable(n :: Node) 
    return Nodes.has_inputs(n.op)
end

function shrink(n :: Node)
    if !shrinkable(n) 
        return [n]
    end

    return n.inputs
end

# FRep tests

frep_tests = Quickcheck("FRep Tests")

function compare(a, b)
    tol = 0.001
    return abs(a - b) <= tol
end

@add_predicate frep_tests "Apply to axes" ((n :: Node, arg :: Vec3{Float64}) -> 
    compare(n(arg), n(Shapes.axes)(arg)))

@add_predicate frep_tests "Constant fold" ((n :: Node, arg :: Vec3{Float64}) -> 
    compare(n(arg), constant_fold(n)(arg)))

@add_predicate frep_tests "Node to tape" ((n :: Node, arg :: Vec3{Float64}) ->
    compare(n(arg), Tapes.run(Tapes.node_to_tape(n), arg.x, arg.y, arg.z)))

@quickcheck frep_tests


