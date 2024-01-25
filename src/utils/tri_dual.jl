# Dual numbers with three partial derivatives.
module TriDuals

export TriDual

struct TriDual <: Number
    val :: Float64
    dx :: Float64
    dy :: Float64
    dz :: Float64
end

# Convert a float constant to a tri-dual.
function Base.convert(::Type{TriDual}, c :: Float64) :: TriDual 
    return TriDual(c, 0., 0., 0.)
end

# This allows for mixed arithmetic between floats and tri-duals.
Base.promote_rule(::Type{TriDual}, ::Type{Float64}) = TriDual

# Basic arithmetic

Base.zero(::Type{TriDual}) = Tridual(0., 0., 0., 0.)

Base.one(::Type{TriDual}) = Tridual(1., 0., 0., 0.)

function Base.:+(a :: TriDual, b :: TriDual) :: TriDual
    return TriDual(
        a.val + b.val,
        a.dx + b.dx,
        a.dy + b.dy,
        a.dz + b.dz)
end

function Base.:-(a :: TriDual) :: TriDual
    return TriDual(
        -a.val,
        -a.dx,
        -a.dy,
        -a.dz)
end

function Base.:-(a :: TriDual, b :: TriDual) :: TriDual
    return TriDual(
        a.val - b.val,
        a.dx - b.dx,
        a.dy - b.dy,
        a.dz - b.dz)
end

function Base.:*(a :: TriDual, b :: TriDual) :: TriDual
    return TriDual(
        a.val * b.val,
        a.dx * b.val + a.val * b.dx,
        a.dy * b.val + a.val * b.dy,
        a.dz * b.val + a.val * b.dz)
end

function Base.inv(a :: TriDual) :: TriDual 
    return TriDual(
        1. / a.val,
        -a.dx / a.val^2,
        -a.dy / a.val^2,
        -a.dz / a.val^2)
end

function Base.:/(a :: TriDual, b :: TriDual) :: TriDual
    return a * inv(b)
end

# More elaborate arithmetic

function Base.sqrt(a :: TriDual) :: TriDual
    return TriDual(
        sqrt(a.val),
        a.dx / (2. * sqrt(a.val)),
        a.dy / (2. * sqrt(a.val)),
        a.dz / (2. * sqrt(a.val)))
end

function Base.sin(a :: TriDual) :: TriDual
    return TriDual(
        sin(a.val),
        a.dx * cos(a.val),
        a.dy * cos(a.val),
        a.dz * cos(a.val))
end

function Base.cos(a :: TriDual) :: TriDual
    return TriDual(
        cos(a.val),
        -a.dx * sin(a.val),
        -a.dy * sin(a.val),
        -a.dz * sin(a.val))
end

function Base.exp(a :: TriDual) :: TriDual
    return TriDual(
        exp(a.val),
        a.dx * exp(a.val),
        a.dy * exp(a.val),
        a.dz * exp(a.val))
end

# Min and max

function Base.min(a :: TriDual, b :: TriDual) :: TriDual
    return a.val <= b.val ? a : b
end

function Base.max(a :: TriDual, b :: TriDual) :: TriDual
    return a.val >= b.val ? a : b
end

end 