

# The log-sum-exp function is a smooth approximation to the maximum.
# >>> log-sum-exp(x, y) = log(exp(x) + exp(y))
# We use a trick to avoid overflows.
function log_sum_exp(x, y)
    if x >= y
        x + log((1. + exp(y - x)) / 1.)
    else 
        y + log((exp(x - y) + 1.) / 1.)
    end
end

# Although this is called the softmax, it is NOT an approximation to the max.
# Rather, it approximates the argmax function, and it shows up as the 
# gradient of log_sum_exp.
# >>> softmax(x, y) = ( exp(x) / (exp(x) + exp(y)) , exp(y) / (exp(x) + exp(y)) )
function softmax(x, y)
    if x >= y
        (1. / (1. + exp(y - x)), exp(y - x) / (1. + exp(y - x)))
    else
        (exp(x - y) / (exp(x - y) + 1.), 1. / (exp(x - y) + 1.))
    end
end

# To obtain a smooth maximum from the log-sum-exp, we use a scaling factor
# to determine how sharp the approximation should be.
const SMAX_SCALE = 20.
@inline smooth_max(x, y) = log_sum_exp(SMAX_SCALE * x, SMAX_SCALE * y) / SMAX_SCALE

# Partial derivatives of the smooth_max function.
@inline smooth_max_grad(x, y) = softmax(SMAX_SCALE * x, SMAX_SCALE * y)

# Same thing as smooth_max, but the scaling factor is negative.
const SMIN_SCALE = -20.
@inline smooth_min(x, y) = log_sum_exp(SMIN_SCALE * x, SMIN_SCALE * y) / SMIN_SCALE

# Partial derivatives of the smooth_min function.
@inline smooth_min_grad(x, y) = softmax(SMIN_SCALE * x, SMIN_SCALE * y)
