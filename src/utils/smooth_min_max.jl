module SmoothMinMax

export log_sum_exp, softmax, smooth_max, smooth_max_grad, smooth_min, smooth_min_grad

# The log-sum-exp function is a smooth approximation to the maximum.
# >>> log-sum-exp(x, y) = log(exp(x) + exp(y))
# We use a trick to avoid overflows.
function log_sum_exp(x :: Float64, y :: Float64)
    if x >= y
        x + log(1. + exp(y - x))
    else 
        y + log(exp(x - y) + 1.)
    end
end

# Although this is called the softmax, it is NOT an approximation to the max.
# Rather, it approximates the argmax function, and it shows up as the 
# gradient of log_sum_exp.
# >>> softmax(x, y) = ( exp(x) / (exp(x) + exp(y)) , exp(y) / (exp(x) + exp(y)) )
function softmax(x :: Float64, y :: Float64)
    if x >= y
        (1. / (1. + exp(y - x)), exp(y - x) / (1. + exp(y - x)))
    else
        (exp(x - y) / (exp(x - y) + 1.), 1. / (exp(x - y) + 1.))
    end
end

# To obtain a smooth maximum from the log-sum-exp, we use a scaling factor
# to determine how sharp the approximation should be.
const SCALE = 1.
@inline function smooth_max(x :: Float64, y :: Float64, scale :: Float64 = SCALE)
    @assert scale > 0.
    log_sum_exp(scale * x, scale * y) / scale
end

# Partial derivatives of the smooth_max function.
@inline function smooth_max_grad(x :: Float64, y :: Float64, scale :: Float64 = SCALE)
    @assert scale > 0.
    softmax(scale * x, scale * y)
end

# Same thing as smooth_max, but the scaling factor is negative.
@inline function smooth_min(x :: Float64, y :: Float64, scale :: Float64 = -SCALE)
    @assert scale < 0.
    log_sum_exp(scale * x, scale * y) / scale
end

# Partial derivatives of the smooth_min function.
@inline function smooth_min_grad(x :: Float64, y :: Float64, scale :: Float64 = -SCALE)
    @assert scale < 0.
    softmax(scale * x, scale * y)
end

end