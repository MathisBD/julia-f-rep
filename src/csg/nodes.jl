module Nodes
using ..Vec3s, ..DotGraph
import Base: convert, show
import ..SmoothMinMax: smooth_max, smooth_min

export Op, Node, X, Y, Z, Const, Neg, Sin, Cos, Exp, Sqrt, Add, Sub, Mul, Div, Min, Max, SMin, SMax
export topoiter, topomap, visualize, constant_fold, merge_axes

@enum Op X Y Z Const Neg Sin Cos Exp Sqrt Add Sub Mul Div Min Max SMin SMax

struct Node
    op :: Op
    inputs :: Vector{Node}
    constant :: Float64
 
    Node(c :: Float64) = new(Const, [], c)
    Node(op :: Op) = new(op, [], 0.0)
    Node(op :: Op, inp1 :: Node) = new(op, [inp1], 0.0)
    Node(op :: Op, inp1 :: Node, inp2 :: Node) = new(op, [inp1, inp2], 0.0)
    #Node(op :: Op, inp1 :: Node, inp2 :: Node, c :: Float64) = new(op, [inp1, inp2], c)
end


function has_one_input(op :: Op) :: Bool
    return op == Neg || op == Sin || op == Cos || op == Exp || op == Sqrt
end

function has_two_inputs(op :: Op) :: Bool
    return op == Add || op == Sub || op == Mul || op == Div || op == Min || op == Max || op == SMin || op == SMax
end

function has_inputs(op :: Op) :: Bool
    return has_one_input(op) || has_two_inputs(op)
end


function Base.show(io :: IO, n :: Node)
    if n.op == X || n.op == Y || n.op == Z
        print(io, "Node($(n.op))")
    elseif n.op == Const
        print(io, "Node($(n.constant))")
    elseif has_one_input(n.op)
        print(io, "Node($(n.op), $(n.inputs[1]))")
    elseif has_two_inputs(n.op) 
        #if !has_const(n.op)
            print(io, "Node($(n.op), $(n.inputs[1]), $(n.inputs[2]))")
        #else 
        #    print(io, "Node($(n.op), $(n.inputs[1]), $(n.inputs[2]), $(n.constant))")
        #end
    else 
        error("Unhandled op $(n.op) :: $(typeof(n.op))")
    end
end

function is_constant(node :: Node, constant :: Float64)
    return node.op == Const && node.constant == constant
end


# Apply an operator to numeric arguments. 
# This doesn't work for all operators.
function apply_op(::Val{Neg}, a :: T) :: T where {T} 
    return -a
end
function apply_op(::Val{Sin}, a :: T) :: T where {T} 
    return sin(a)
end
function apply_op(::Val{Cos}, a :: T) :: T where {T} 
    return cos(a)
end
function apply_op(::Val{Exp}, a :: T) :: T where {T} 
    return exp(a)
end
function apply_op(::Val{Sqrt}, a :: T) :: T where {T} 
    return sqrt(a)
end
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
function apply_op(::Val{SMin}, a :: T, b :: T) :: T where {T} 
    return smooth_min(a, b)
end
function apply_op(::Val{SMax}, a :: T, b :: T) :: T where {T} 
    return smooth_max(a, b)
end

# Convert a float to a node.
convert(::Type{Node}, x :: Float64) = Node(x)

Base.:-(a :: Node) = Node(Neg, a)
Base.sin(a :: Node) = Node(Sin, a)
Base.cos(a :: Node) = Node(Cos, a)
Base.exp(a :: Node) = Node(Exp, a)
Base.sqrt(a :: Node) = Node(Sqrt, a)

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


smooth_min(a :: Node, b :: Node)    = Node(SMin, a, b)
smooth_min(a :: Node, b :: Float64) = Node(SMin, a, Node(b))
smooth_min(a :: Float64, b :: Node) = Node(SMin, Node(a), b)

smooth_max(a :: Node, b :: Node)    = Node(SMax, a, b)
smooth_max(a :: Node, b :: Float64) = Node(SMax, a, Node(b))
smooth_max(a :: Float64, b :: Node) = Node(SMax, Node(a), b)

#smooth_min(a :: Node, b :: Node, scale :: Float64)    = Node(SMin, a, b, scale)
#smooth_min(a :: Node, b :: Float64, scale :: Float64) = Node(SMin, a, Node(b), scale)
#smooth_min(a :: Float64, b :: Node, scale :: Float64) = Node(SMin, Node(a), b, scale)
#
#smooth_max(a :: Node, b :: Node, scale :: Float64)    = Node(SMax, a, b, scale)
#smooth_max(a :: Node, b :: Float64, scale :: Float64) = Node(SMax, a, Node(b), scale)
#smooth_max(a :: Float64, b :: Node, scale :: Float64) = Node(SMax, Node(a), b, scale)

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
        #elseif n.op == SMin
        #    return smooth_min(children..., n.constant)
        #elseif n.op == SMax
        #    return smooth_max(children..., n.constant)
        else
            return apply_op(Val(n.op), children...)
        end
    end

    return topomap(T, helper, root)
end

function constant_fold_step(node :: Node, new_inputs :: Vector{Node})
    # The node has no inputs : we can't do anything.
    if !has_inputs(node.op)
        return node
    end

    # All the new inputs are constant : return a constant node.
    constants = [inp.constant for inp in new_inputs if inp.op == Const]
    if length(constants) == length(new_inputs)
        return Node(apply_op(Val(node.op), constants...))
    end

    # Some inputs are constants : we can apply some tricks.
    if node.op == Add
        if is_constant(new_inputs[1], 0.0) 
            return new_inputs[2]
        elseif is_constant(new_inputs[2], 0.0) 
            return new_inputs[1]
        end
    elseif node.op == Sub
        if is_constant(new_inputs[2], 0.0)
            return new_inputs[1]
        end
    elseif node.op == Mul
        if is_constant(new_inputs[1], 0.0) || is_constant(new_inputs[2], 0.0)
            return Node(0.0)
        elseif is_constant(new_inputs[1], 1.0)
            return new_inputs[2]
        elseif is_constant(new_inputs[2], 1.0)
            return new_inputs[1]
        end
    elseif node.op == Div
        if is_constant(new_inputs[1], 0.0) 
            return Node(0.0)
        elseif is_constant(new_inputs[2], 1.0)
            return new_inputs[1]
        elseif new_inputs[2].op == Const
            return new_inputs[1] * (1.0 / new_inputs[2].constant)
        end
    end

    # We ran out of tricks : return the node with the new inputs.
    return Node(node.op, new_inputs...)
end


function constant_fold(root :: Node)
    return topomap(Node, constant_fold_step, root)
end


function merge_axes(root :: Node)
    local x :: Union{Nothing, Node} = nothing
    local y :: Union{Nothing, Node} = nothing
    local z :: Union{Nothing, Node} = nothing
    
    function helper(node :: Node, new_inputs :: Vector{Node})
        if node.op == X
            if isnothing(x) 
                x = node
            end
            return x
        elseif node.op == Y
            if isnothing(y) 
                y = node
            end
            return y
        elseif node.op == Z
            if isnothing(z) 
                z = node
            end
            return z
        elseif node.op == Const
            return Node(node.constant)
        else 
            return Node(node.op, new_inputs...)
        end
    end

    return topomap(Node, helper, root)
end

end


