# This is inside module Nodes.

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