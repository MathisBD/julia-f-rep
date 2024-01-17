# This is inside module Nodes.

function constant_fold_step(node :: Node{T}, new_inputs :: Vector{Node{T}}) where {T}
    # The node has no inputs : we can't do anything.
    if !has_inputs(node.op)
        return node
    end

    # All the new inputs are constant : return a constant node.
    constants = [inp.constant for inp in new_inputs if inp.op == Const]
    if length(constants) == length(new_inputs)
        return Node{T}(apply_op(Val(node.op), constants...))
    end

    # Some inputs are constants : we can apply some tricks.
    if node.op == Add
        if is_constant(new_inputs[1], zero(T)) 
            return new_inputs[2]
        elseif is_constant(new_inputs[2], zero(T)) 
            return new_inputs[1]
        end
    elseif node.op == Sub
        if is_constant(new_inputs[2], zero(T))
            return new_inputs[1]
        end
    elseif node.op == Mul
        if is_constant(new_inputs[1], zero(T)) || is_constant(new_inputs[2], zero(T))
            return Node{T}(zero(T))
        elseif is_constant(new_inputs[1], one(T))
            return new_inputs[2]
        elseif is_constant(new_inputs[2], one(T))
            return new_inputs[1]
        end
    elseif node.op == Div
        if is_constant(new_inputs[1], zero(T)) 
            return Node{T}(zero(T))
        elseif is_constant(new_inputs[2], one(T))
            return new_inputs[1]
        elseif new_inputs[2].op == Const
            return new_inputs[1] * (one(T) / new_inputs[2].constant)
        end
    end

    # We ran out of tricks : return the node with the new inputs.
    return Node{T}(node.op, new_inputs...)
end


function constant_fold(root :: Node{T}) where {T}
    return topomap(Node{T}, constant_fold_step, root)
end