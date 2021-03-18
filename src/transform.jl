function translate_ray(ray::Ray{T}, xoff::T, yoff::T, zoff::T) where {T<:AbstractFloat}
    x = ray.x + xoff
    y = ray.y + yoff
    z = ray.z + zoff
    tray = Ray(x, y, z, ray.dx, ray.dy, ray.dz)
    return tray
end

function transform_tzy(ray::Ray{T}, xoff::T, yoff::T, zoff::T, beta::T, gamma::T) where {T<:AbstractFloat}
    # to transform from world space into object space
    # (values for xoff, yoff, ..., gamma, should be of opposite sign to how they are defined in the scene specification file)
    
    # the translation
    x = ray.x + xoff
    y = ray.y + yoff
    z = ray.z + zoff
    
    # rotate about z-axis
    x_ = x * cos(gamma) + y * -sin(gamma)
    y  = x * sin(gamma) + y * cos(gamma)
    dx_ = ray.dx * cos(gamma) + ray.dy * -sin(gamma)
    dy  = ray.dx * sin(gamma) + ray.dy * cos(gamma)
    
    # rotate about y-axis
    x = x_ * cos(beta) + z * sin(beta)
    z = x_ * -sin(beta) + z * cos(beta)
    dx = dx_ * cos(beta) + ray.dz * sin(beta)
    dz = dx_ * -sin(beta) + ray.dz * cos(beta)
    
    tray = Ray(x, y, z, dx, dy, dz)
    return tray
end

function carthesian2spherical(x::T, y::T, z::T) where {T<:AbstractFloat}
    theta = atan(sqrt(x*x + y*y), z)
    phi = atan(x, y)
    if phi < 0.0
        phi += 2.0*pi
    end
    return theta, phi
end

function spherical2carthesian(theta::T, phi::T) where {T<:AbstractFloat}
    x = sin(theta) * sin(phi)
    y = sin(theta) * cos(phi)
    z = cos(theta)
    return x, y, z
end

function relazi(vx::T, vy::T, vz::T, nx::T, ny::T, nz::T, sx::T, sy::T, sz::T) where {T<:AbstractFloat}
    # v=viewVector; n=normalVector; h=halfVector...
    ux, uy, uz = cross_scalars(nx, ny, nz, sx, sy, sz)
    vx, vy, vz = cross_scalars(nx, ny, nz, vx, vy, vz)
    nsx, nsy, nsz = normalize_scalars(ux, uy, uz)
    nvx, nvy, nvz = normalize_scalars(vx, vy, vz)
    a = acos( nsx*nvx + nsy*nvy + nsz*nvz )
    return a
end

function get_extremes(x::Vector{T}, y::Vector{T}, z::Vector{T}) where {T<:AbstractFloat}
    minx, maxx = minimum(x), maximum(x)
    miny, maxy = minimum(y), maximum(y)
    minz, maxz = minimum(z), maximum(z)
    p::Vector{T} = [minx, minx, minx, minx, maxx, maxx, maxx, maxx]
    q::Vector{T} = [miny, miny, maxy, maxy, miny, miny, maxx, maxx]
    r::Vector{T} = [minz, maxz, minz, maxz, minz, maxz, minz, maxz]
    return p, q, r
end

function get_extremes(geometry::Geometry{T,S,Array{T,1},Array{S,1}}) where {T<:AbstractFloat, S<:Integer}
    x = vcat(geometry.v1x, geometry.v2x, geometry.v3x)
    y = vcat(geometry.v1y, geometry.v2y, geometry.v3y)
    z = vcat(geometry.v1z, geometry.v2z, geometry.v3z)
    p, q, r = get_extremes(x, y, z)
    return p, q, r
end

function get_extremes(geometry::Geometry{T,S,A,E}, hook::Hook{T,S}) where {T<:AbstractFloat, S<:Integer, A<:AbstractArray{T,1}, E<:AbstractArray{S,1}}
    x, y, z = get_extremes(geometry)
    x, y, z = transform_forward(x, y, z, hook)
    x, y, z = get_extremes(x, y, z)
    return x, y, z
end

function transform(x::T, y::T, z::T, xoff::T, yoff::T, zoff::T, beta::T, gamma::T) where {T<:AbstractFloat}
    # to transform x,y,z from object space into worldspace
    # rotate about y-axis
    t = x * cos(beta) + z * sin(beta)
    c = x * -sin(beta) + z * cos(beta)

    # rotate about z-axis
    a = t * cos(gamma) + y * -sin(gamma)
    b = t * sin(gamma) + y * cos(gamma)

    # and the translation
    a = a + xoff
    b = b + yoff
    c = c + zoff
    return a, b, c
end

function translate_forward(x::T, y::T, z::T, hook::Hook{T,S}) where {T<:AbstractFloat, S<:Integer}
    x = x + hook.xoff
    y = y + hook.yoff
    z = z + hook.zoff
    return x, y, z
end

function translate_inverse(x::T, y::T, z::T, hook::Hook{T,S}) where {T<:AbstractFloat, S<:Integer}
    x = x - hook.xoff
    y = y - hook.yoff
    z = z - hook.zoff
    return x, y, z
end

function translate_forward(x::AbstractArray{T}, y::AbstractArray{T}, z::AbstractArray{T}, hook::Hook{T,S}) where {T<:AbstractFloat, S<:Integer}
    x = x .+ hook.xoff
    y = y .+ hook.yoff
    z = z .+ hook.zoff
    return x, y, z
end

function translate_inverse(x::AbstractArray{T}, y::AbstractArray{T}, z::AbstractArray{T}, hook::Hook{T,S}) where {T<:AbstractFloat, S<:Integer}
    x = x .- hook.xoff
    y = y .- hook.yoff
    z = z .- hook.zoff
    return x, y, z
end

function rotate_forward(x::T, y::T, z::T, hook::Hook{T,S}) where {T <: AbstractFloat, S<:Integer}
    # rotate about y-axis
    nx_ = x * cos(hook.beta) + z * sin(hook.beta)
    nz =  x * -1.0 * sin(hook.beta) + z * cos(hook.beta)

    # rotate about z-axis
    nx =  nx_ * cos(hook.gamma) + y * -1.0 * sin(hook.gamma)
    ny =  nx_ * sin(hook.gamma) + y * cos(hook.gamma)
    return nx, ny, nz
end

function rotate_forward(x::AbstractArray{T}, y::AbstractArray{T}, z::AbstractArray{T}, hook::Hook{T,S}) where {T <: AbstractFloat, S<:Integer}
    # rotate about y-axis
    nx_ = x .* cos(hook.beta) .+ z .* sin(hook.beta)
    nz =  x .* -1.0 .* sin(hook.beta) .+ z .* cos(hook.beta)

    # rotate about z-axis
    nx =  nx_ .* cos(hook.gamma) .+ y .* -1.0 .* sin(hook.gamma)
    ny =  nx_ .* sin(hook.gamma) .+ y .* cos(hook.gamma)
    return nx, ny, nz
end

function rotate_inverse(x::T, y::T, z::T, hook::Hook{T,S}) where {T<:AbstractFloat, S<:Integer}
    inv_beta = -1.0 * hook.beta
    inv_gamma = -1.0 * hook.gamma

    # rotate about z-axis (to get segment tip on y=0 plane)
    nx_1 = x * cos(inv_gamma) + y * -1.0 * sin(inv_gamma)
    ny = x * sin(inv_gamma) + y * cos(inv_gamma)

    # rotate about y-axis (to get segment tip on z=0 plane)
    nx = nx_1 * cos(inv_beta) + z * sin(inv_beta)
    nz = nx_1 * -1.0 * sin(inv_beta) + z * cos(inv_beta)
    return nx, ny, nz
end

function transform_forward(x::AbstractArray{T}, y::AbstractArray{T}, z::AbstractArray{T}, hook::Hook{T,S}) where {T<:AbstractFloat, S<:Integer}
    x, y, z = rotate_forward(x, y, z, hook)
    x, y, z = translate_forward(x, y, z, hook)
    return x, y, z
end

function transform_inverse(ray::Ray{T}, hook::Hook{T,S}) where {T<:AbstractFloat, S<:Integer}
    x, y, z = translate_inverse(ray.x, ray.y, ray.z, hook)
    x, y, z = rotate_inverse(x, y, z, hook)
    dx, dy, dz = rotate_inverse(ray.dx, ray.dy, ray.dz, hook)
    transRay::Ray = Ray(x, y, z, dx, dy, dz)
    return transRay
end


function weibull(x, a, b)
    y = (a/b) .* (x./b).^(a-1) .* exp.(.-(x./b).^a)
    ys  = y ./ maximum(y)
    return y, ys
end

