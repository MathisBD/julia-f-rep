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
        if is_constant(new_inputs[1], 0) 
            return new_inputs[2]
        elseif is_constant(new_inputs[2], 0) 
            return new_inputs[1]
        end
    elseif node.op == Sub
        if is_constant(new_inputs[2], 0)
            return new_inputs[1]
        end
    elseif node.op == Mul
        if is_constant(new_inputs[1], 0) || is_constant(new_inputs[2], 0)
            return Node{T}(0)
        elseif is_constant(new_inputs[1], 1)
            return new_inputs[2]
        elseif is_constant(new_inputs[2], 1)
            return new_inputs[1]
        end
    end

    # We ran out of tricks : return the node with the new inputs.
    return Node{T}(node.op, new_inputs...)
end


function constant_fold(root :: Node{T}) where {T}
    return topomap(Node{T}, constant_fold_step, root)
end