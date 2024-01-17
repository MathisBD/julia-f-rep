using JCheck
import JCheck: generate, shrinkable, shrink

include("../src/frep.jl")
using .Vec3s, .Nodes, .Shapes

# FRep generators

#function generate(rng, ::Type{Op}, n :: Int)
#    return rand(rng, instances(Op), n)
#end

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
        #op = rand(rng, (Add, Sub, Mul, Div, Min, Max))
        op = rand(rng, (Add, Sub, Mul, Min, Max))
    end

    # Axis node.
    if op == X || op == Y || op == Z 
        @assert size == 0
        return Node{Float64}(op)
    # Constant node.
    elseif op == Const
        @assert size == 0
        return Node{Float64}(rand_float64(rng))
    # Node with two inputs
    else
        @assert size > 0
        
        size_left = rand(rng, 0:(size-1))
        size_right = size - 1 - size_left
        
        left = rand_node_sized(rng, size_left)
        right = rand_node_sized(rng, size_right)
        return Node{Float64}(op, left, right) 
    end
end

function generate(rng, ::Type{Node{Float64}}, n :: Int)
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

function shrinkable(n :: Node{Float64}) 
    return Nodes.has_inputs(n.op)
end

function shrink(n :: Node{Float64})
    if !shrinkable(n) 
        return [n]
    end

    return n.inputs
end

# FRep tests

frep_tests = Quickcheck("FRep Tests")

function compare(a :: Node{T}, b :: Node{T}, arg :: Vec3{T}) where {T}
    tol = 0.001
    return abs(a(arg) - b(arg)) <= tol
end

@add_predicate frep_tests "Apply to axes" ((n :: Node{Float64}, arg :: Vec3{Float64}) -> 
    compare(n, n(Shapes.axes(Float64)), arg))

@add_predicate frep_tests "Constant fold" ((n :: Node{Float64}, arg :: Vec3{Float64}) -> 
    compare(n, constant_fold(n), arg))

@quickcheck frep_tests


