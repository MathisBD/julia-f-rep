module Intervals

export Interval

# The bounds of the interval should not be NaNs,
# but they can be infinite (-Inf for the lower bound and Inf for the upper bound).
# Both bounds can not be equal to the same infinity.
struct Interval <: Number
    low :: Float64
    high :: Float64
end

function is_valid(a :: Interval) :: Bool
    return !isnan(a.low) && !isnan(a.high) && a.low != Inf && a.high != -Inf && a.low <= a.high 
end

# Handles NaNs in the bounds. 
function safe_interval(low :: Float64, high :: Float64) :: Interval
    return Interval(
        isnan(low) ? -Inf : low,
        isnan(high) ? Inf : high)
end

# Convert a float constant to an interval.
function Base.convert(::Type{Interval}, c :: Float64) :: Interval 
    return safe_interval(c, c)
end

# This allows for mixed arithmetic between floats and intervals.
Base.promote_rule(::Type{Interval}, ::Type{Float64}) = Interval

# Basic arithmetic

function Base.:-(a :: Interval) :: Interval
    return Interval(
        -a.high,
        -a.low)
end

function Base.:+(a :: Interval, b :: Interval) :: Interval
    return Interval(
        a.low + b.low,
        a.high + b.high)
end

function Base.:-(a :: Interval, b :: Interval) :: Interval
    return Interval(
        a.low - b.high,
        a.high - b.low)
end

function Base.:*(a :: Interval, b :: Interval) :: Interval
    return safe_interval(
        min(a.low * b.low, a.low * b.high, a.high * b.low, a.high * b.high),
        max(a.low * b.low, a.low * b.high, a.high * b.low, a.high * b.high))
end

function Base.inv(a :: Interval) :: Interval
    if a.low <= 0. <= a.high
        return Interval(-Inf, Inf)
    else
        return safe_interval(1. / a.high, 1. / a.low)
    end
end

function Base.:/(a :: Interval, b :: Interval) :: Interval
    return a * inv(b)
end

# More elaborate arithmetic

function Base.min(a :: Interval, b :: Interval) :: Interval
    return Interval(
        min(a.low, b.low),
        min(a.high, b.high))
end

function Base.max(a :: Interval, b :: Interval) :: Interval
    return Interval(
        max(a.low, b.low),
        max(a.high, b.high))
end

function Base.exp(a :: Interval) :: Interval
    return Interval(
        exp(a.low),
        exp(a.high))
end

function Base.sqrt(a :: Interval) :: Interval
    return Interval(
        0. <= a.low ? sqrt(a.low) : 0., 
        0. <= a.high ? sqrt(a.high) : 0.)        
end

@inline function contains_int(low :: Float64, high :: Float64) :: Bool
    return low <= floor(high)
end

function Base.sin(a :: Interval) :: Interval
    if a.low == -Inf || a.high == Inf
        return Interval(-1., 1.)
    else
        TWO_PI = 2 * pi
        low = contains_int(a.low/TWO_PI - 3/4, a.high/TWO_PI - 3/4) ? 
              -1. : min(sin(a.low), sin(a.high))    
        high = contains_int(a.low/TWO_PI - 1/4, a.high/TWO_PI - 1/4) ?
               1. : max(sin(a.low), sin(a.high))
        return Interval(low, high)
    end
end

function Base.cos(a :: Interval) :: Interval
    if a.low == -Inf || a.high == Inf
        return Interval(-1., 1.)
    else
        TWO_PI = 2 * pi
        low = contains_int(a.low/TWO_PI - 1/2, a.high/TWO_PI - 1/2) ? 
              -1. : min(cos(a.low), cos(a.high))    
        high = contains_int(a.low/TWO_PI, a.high/TWO_PI) ?
               1. : max(cos(a.low), cos(a.high))
        return Interval(low, high)
    end
end

end