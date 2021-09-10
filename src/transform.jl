function weibull(x, a, b)
    y = (a/b) .* (x./b).^(a-1) .* exp.(.-(x./b).^a)
    ys  = y ./ maximum(y)
    return y, ys
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
    dotprod = nsx*nvx + nsy*nvy + nsz*nvz
    if dotprod > 1.0 # it would only be slightly higher by machine precision, e.g. 1.0000000000000002
        dotprod = 1.0
    end
    if dotprod < 0.0 # similar to the case above...
        dotprod = 0.0
    end
    a = acos(dotprod)
    return a
end

function translate(ray::Ray{T}, xoff::T, yoff::T, zoff::T) where {T<:AbstractFloat}
    x = ray.x + xoff
    y = ray.y + yoff
    z = ray.z + zoff
    tray = Ray(x, y, z, ray.dx, ray.dy, ray.dz)
    return tray
end

function transform(ray::Ray{T}, xoff::T, yoff::T, zoff::T, beta::T, gamma::T) where {T<:AbstractFloat}
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



# functions using Hook

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

function rotate_forward(x::AbstractArray{T}, y::AbstractArray{T}, z::AbstractArray{T}, hook::Hook{T,S}) where {T<:AbstractFloat, S<:Integer}
    # rotate about y-axis
    nx_ = x .* cos(hook.beta) .+ z .* sin(hook.beta)
    nz =  x .* -1.0 .* sin(hook.beta) .+ z .* cos(hook.beta)

    # rotate about z-axis
    nx =  nx_ .* cos(hook.gamma) .+ y .* -1.0 .* sin(hook.gamma)
    ny =  nx_ .* sin(hook.gamma) .+ y .* cos(hook.gamma)
    return nx, ny, nz
end

function rotate_inverse(x::AbstractArray{T}, y::AbstractArray{T}, z::AbstractArray{T}, hook::Hook{T,S}) where {T<:AbstractFloat, S<:Integer}
    inv_beta = -1.0 * hook.beta
    inv_gamma = -1.0 * hook.gamma

    # rotate about z-axis (to get segment tip on y=0 plane)
    nx_1 = x .* cos(inv_gamma) .+ y .* -1.0 .* sin(inv_gamma)
    ny = x .* sin(inv_gamma) .+ y .* cos(inv_gamma)

    # rotate about y-axis (to get segment tip on z=0 plane)
    nx = nx_1 .* cos(inv_beta) .+ z .* sin(inv_beta)
    nz = nx_1 .* -1.0 .* sin(inv_beta) .+ z .* cos(inv_beta)
    return nx, ny, nz
end

function transform_forward(x::T, y::T, z::T, hook::Hook{T,S}) where {T<:AbstractFloat, S<:Integer}
    x, y, z = rotate_forward(x, y, z, hook)
    x, y, z = translate_forward(x, y, z, hook)
    return x, y, z
end

function transform_inverse(x::T, y::T, z::T, hook::Hook{T,S}) where {T<:AbstractFloat, S<:Integer}
    x, y, z = translate_inverse(x, y, z, hook)
    x, y, z = rotate_inverse(x, y, z, hook)
    return x, y, z
end

function transform_forward(x::AbstractArray{T}, y::AbstractArray{T}, z::AbstractArray{T}, hook::Hook{T,S}) where {T<:AbstractFloat, S<:Integer}
    x, y, z = rotate_forward(x, y, z, hook)
    x, y, z = translate_forward(x, y, z, hook)
    return x, y, z
end

function transform_inverse(x::AbstractArray{T}, y::AbstractArray{T}, z::AbstractArray{T}, hook::Hook{T,S}) where {T<:AbstractFloat, S<:Integer}
    x, y, z = translate_inverse(x, y, z, hook)
    x, y, z = rotate_inverse(x, y, z, hook)
    return x, y, z
end

function transform_forward(ray::Ray{T}, hook::Hook{T,S}) where {T<:AbstractFloat, S<:Integer}
    x, y, z = rotate_forward(ray.x, ray.y, ray.z, hook)
    x, y, z = translate_forward(x, y, z, hook)
    dx, dy, dz = rotate_forward(ray.dx, ray.dy, ray.dz, hook)
    transRay::Ray = Ray(x, y, z, dx, dy, dz)
    return transRay
end

function transform_inverse(ray::Ray{T}, hook::Hook{T,S}) where {T<:AbstractFloat, S<:Integer}
    x, y, z = translate_inverse(ray.x, ray.y, ray.z, hook)
    x, y, z = rotate_inverse(x, y, z, hook)
    dx, dy, dz = rotate_inverse(ray.dx, ray.dy, ray.dz, hook)
    transRay::Ray = Ray(x, y, z, dx, dy, dz)
    return transRay
end



# functions using transformation matrix

function transform(x::T, y::T, z::T, matrix::AbstractArray{T,2}) where {T<:AbstractFloat}
    xnew = x*matrix[1,1] + y*matrix[1,2] + z*matrix[1,3] + matrix[1,4]
    ynew = x*matrix[2,1] + y*matrix[2,2] + z*matrix[2,3] + matrix[2,4]
    znew = x*matrix[3,1] + y*matrix[3,2] + z*matrix[3,3] + matrix[3,4]
    return xnew, ynew, znew
end

function transform(ray::Ray{T}, matrix::AbstractArray{T,2}) where {T<:AbstractFloat, S<:Integer}
    x2 = ray.x + ray.dx
    y2 = ray.y + ray.dy
    z2 = ray.z + ray.dz
    x, y, z = transform(ray.x, ray.y, ray.z, matrix)
    u, v, w = transform(x2, y2, z2, matrix)
    dx = u - x
    dy = v - y
    dz = w - z
    #dx, dy, dz = normalize(dx, dy, dz)   # DON'T APPLY NORMALIZE !!!
    transformedRay::Ray = Ray(x, y, z, dx, dy, dz)
    return transformedRay
end



# function using Transform (<:Asset) object

function transform_forward(x::T, y::T, z::T, t::Transform{T,S}) where {T<:AbstractFloat, S<:Integer}
    return transform(x, y, z, t.matrix)
end

function transform_inverse(x::T, y::T, z::T, t::Transform{T,S}) where {T<:AbstractFloat, S<:Integer}
    return transform(x, y, z, t.inverse)
end

function transform_forward(ray::Ray{T}, t::Transform{T,S}) where {T<:AbstractFloat, S<:Integer}
    return transform(ray, t.matrix)
end

function transform_inverse(ray::Ray{T}, t::Transform{T,S}) where {T<:AbstractFloat, S<:Integer}
    return transform(ray, t.inverse)
end

