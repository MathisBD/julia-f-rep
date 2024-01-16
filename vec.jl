module Vec3s

export Vec3, dist2, dist

# A vector/point in 3-dimensional space.
# The coordinates do not have to be floats.
struct Vec3{T} 
    x :: T
    y :: T
    z :: T
end

function dist2(a :: Vec3{T}, b :: Vec3{T}) where {T}
    return (a.x - b.x)^2 + (a.y - b.y)^2 + (a.z - b.z)^2
end

function dist(a :: Vec3{T}, b :: Vec3{T}) where {T}
    return sqrt(dist2(a, b))
end

# Basic arithmetic on vectors.
Base.:-(a :: Vec3{T}) where {T} = 
    Vec3(-a.x, -a.y, -a.z)
Base.:+(a :: Vec3{T}, b :: Vec3{T}) where {T} = 
    Vec3(a.x + b.x, a.y + b.y, a.z + b.z)
Base.:-(a :: Vec3{T}, b :: Vec3{T}) where {T} = 
    Vec3(a.x - b.x, a.y - b.y, a.z - b.z)
Base.:*(a :: Vec3{T}, scale :: T) where {T} = 
    Vec3(a.x * scale, a.y * scale, a.z * scale)
Base.:*(scale :: T, a :: Vec3{T}) where {T} = 
    Vec3(scale * a.x, scale * a.y, scale * a.z)
Base.:/(a :: Vec3{T}, scale :: T) where {T} = 
    Vec3(a.x / scale, a.y / scale, a.z / scale)

end
