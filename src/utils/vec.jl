module Vec3s

export Vec3, dist2, dist, norm2, norm, normalize, fmap

# A vector/point in 3-dimensional space.
# The coordinates do not have to be floats.
struct Vec3{T} 
    x :: T
    y :: T
    z :: T
end

# If we can convert vector fields, we can convert vectors.
@inline function Base.convert(:: Type{Vec3{T}}, v :: Vec3{U}) where {T, U} 
    return Vec3{T}(convert(T, v.x), convert(T, v.y), convert(T, v.z))
end

@inline function full(a :: T) :: Vec3{T} where {T}
    return Vec3(a, a, a)
end

@inline function fmap(f :: Function, a :: Vec3{T}) where {T}
    return Vec3(f(a.x), f(a.y), f(a.z))
end    

@inline function fmap(X :: Type, a :: Vec3{T}) where {T}
    return Vec3(X(a.x), X(a.y), X(a.z))
end     

@inline function fmap(f :: Function, a :: Vec3{T}, b :: Vec3{T}) where {T}
    return Vec3(f(a.x, b.x), f(a.y, b.y), f(a.z, b.z))
end


# Basic arithmetic on vectors. Operations are done element-wise.
Base.:-(a :: Vec3{T}) where {T} = Vec3(-a.x, -a.y, -a.z)
Base.:+(a :: Vec3{T}, b :: Vec3{T}) where {T} = Vec3(a.x + b.x, a.y + b.y, a.z + b.z)
Base.:-(a :: Vec3{T}, b :: Vec3{T}) where {T} = Vec3(a.x - b.x, a.y - b.y, a.z - b.z)
Base.:*(a :: Vec3{T}, b :: Vec3{T}) where {T} = Vec3(a.x * b.x, a.y * b.y, a.z * b.z)
Base.:/(a :: Vec3{T}, b :: Vec3{T}) where {T} = Vec3(a.x / b.x, a.y / b.y, a.z / b.z)

# Mixed arithmetic (vector and scalar).
Base.:+(a :: Vec3{T}, b :: T) where {T} = a + full(b)
Base.:+(a :: T, b :: Vec3{T}) where {T} = full(a) + b
Base.:-(a :: Vec3{T}, b :: T) where {T} = a - full(b)
Base.:-(a :: T, b :: Vec3{T}) where {T} = full(a) - b
Base.:*(a :: Vec3{T}, b :: T) where {T} = a * full(b)
Base.:*(a :: T, b :: Vec3{T}) where {T} = full(a) * b
Base.:/(a :: Vec3{T}, b :: T) where {T} = a / full(b)
Base.:/(a :: T, b :: Vec3{T}) where {T} = full(a) / b
 
# Comparisons on vectors, element-wise.
Base.:(<=)(a :: Vec3{T}, b :: Vec3{T}) where {T} = Vec3(a.x <= b.x, a.y <= b.y, a.z <= b.z)
Base.:<(a :: Vec3{T}, b :: Vec3{T})    where {T} = Vec3(a.x < b.x, a.y < b.y, a.z < b.z)
@inline Base.:(==)(a :: Vec3{T}, b :: Vec3{T}) where {T} = Vec3(a.x == b.x, a.y == b.y, a.z == b.z)

# Mixed comparisons (vector and scalar).
Base.:(<=)(a :: Vec3{T}, b :: T) where {T} = a <= full(b)
Base.:<(a :: Vec3{T}, b :: T)    where {T} = a < full(b)
@inline Base.:(==)(a :: Vec3{T}, b :: T) where {T} = a == full(b)

Base.:(<=)(a :: T, b :: Vec3{T}) where {T} = full(a) <= b
Base.:<(a :: T, b :: Vec3{T})    where {T} = full(a) < b
Base.:(==)(a :: T, b :: Vec3{T}) where {T} = full(a) == b

# A conditional expression coordinate-wise.
@inline function ite(b :: Vec3{Bool}, v1 :: Vec3{T}, v2 :: Vec3{T}) where {T}
    return Vec3(
        b.x ? v1.x : v2.x, 
        b.y ? v1.y : v2.y, 
        b.z ? v1.z : v2.z)
end

@inline function Base.all(a :: Vec3{Bool})
    return all((a.x, a.y, a.z))
end

@inline function Base.any(a :: Vec3{Bool})
    return any((a.x, a.y, a.z))
end

@inline function Base.minimum(a :: Vec3{T}) :: T where {T}
    return min(a.x, a.y, a.z)
end

@inline function Base.maximum(a :: Vec3{T}) :: T where {T}
    return max(a.x, a.y, a.z)
end

# Absolute value coordiante-wise.
@inline function Base.abs(a :: Vec3{T}) where {T}
    return Vec3(abs(a.x), abs(a.y), abs(a.z))
end

@inline function normalize(a :: Vec3{T}) where {T}
    return a / norm(a)
end

@inline function dist2(a :: Vec3{T}, b :: Vec3{T}) where {T}
    return (a.x - b.x)^2 + (a.y - b.y)^2 + (a.z - b.z)^2
end

@inline function dist(a :: Vec3{T}, b :: Vec3{T}) where {T}
    return sqrt((a.x - b.x)^2 + (a.y - b.y)^2 + (a.z - b.z)^2)
end

@inline function norm2(a :: Vec3{T}) where {T}
    return a.x^2 + a.y^2 + a.z^2
end

@inline function norm(a :: Vec3{T}) where {T}
    return sqrt(a.x^2 + a.y^2 + a.z^2)
end

end
