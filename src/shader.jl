#########################################################################################################
#
#  The following shader functions have inputs that are all in the form (vx, vy, vz, nx, ny, nz, args...)
#
#########################################################################################################




#######################################################
#
# Photon shaders... (i.e. single waveband)
#
#######################################################

function lambertian_shader(Ix::T, Iy::T, Iz::T, vx::T, vy::T, vz::T, nx::T, ny::T, nz::T ; rho::T) where {T<:AbstractFloat}
    # note that vx, vy, and vz are not used but only exist so that input arguments are of the same pattern as in other shader functions...
    t, p = sample_f() # could be replaced with a function RPVSampler if we knew how to implement importance sampling for RPV...
    fx, fy, fz = spherical2carthesian(t, p)
    sx, sy, sz = project2normal(fx, fy, fz, nx, ny, nz)
    newray = propagate(Ix, Iy, Iz, sx, sy, sz)
    return true, newray, rho
end

function bilambertian_shader(Ix::T, Iy::T, Iz::T, vx::T, vy::T, vz::T, nx::T, ny::T, nz::T ; rho::T, tau::T) where {T<:AbstractFloat}
    # note that vx, vy, and vz are not used but only exist so that input arguments are of the same pattern as in other shader functions...
    l = rho + tau
    r::Float64
    if rand() > (rho/l)
        nx, ny, nz = -1.0*nx, -1.0*ny, -1.0*nz
        r = tau
    else
        r = rho
    end
    t, p = sample_f() # could be replaced with a function RPVSampler if we knew how to implement importance sampling for RPV...
    fx, fy, fz = spherical2carthesian(t, p)
    sx, sy, sz = project2normal(fx, fy, fz, nx, ny, nz)
    newray = propagate(Ix, Iy, Iz, sx, sy, sz)
    return true, newray, r
end

function pure_reflection_shader(vx::T, vy::T, vz::T, nx::T, ny::T, nz::T ; rho::T) where {T<:AbstractFloat}
    sx, sy, sz = pure_reflection(vx, vy, vz, nx, ny, nz)
    newray = propagate(Ix, Iy, Iz, sx, sy, sz)
    return true, newray, rho
end

function rpv_shader(vx::T, vy::T, vz::T, nx::T, ny::T, nz::T ; rho_zero::T, rho_c::T, captheta::T, k::T) where {T<:AbstractFloat}
    t, p = sample_f() # could be replaced with a function RPVSampler if we knew how to implement importance sampling for RPV...
    fx, fy, fz = spherical2carthesian(t, p)
    sx, sy, sz = project2world(fx, fy, fz, nx, ny, nz)
    phi = relazi(vn, vy, vz, nx, ny, nz, sx, sy, sz)
    theta = acos(vx*nx + vy*ny + vz*nz)
    theta_zero = acos(sx*nx + sy*ny + sz*nz)
    rho_sfc = rpv(theta_zero, theta, phi, rho_zero, rho_c, captheta, k)
    newray = propagate(Ix, Iy, Iz, sx, sy, sz)
    return true, newray, rho_sfc
end




###############################################################
#
#  Path shaders... (i.e. acting on spectra)
#
###############################################################

function lambertian_path_shader(Ix::T, Iy::T, Iz::T, vx::T, vy::T, vz::T, nx::T, ny::T, nz::T ; rho::Array{T}) where {T<:AbstractFloat}
    # note that vx, vy, and vz are not used but only exist so that input arguments are of the same pattern as in other shader functions...
    t, p = sample_f() # could be replaced with a function RPVSampler if we knew how to implement importance sampling for RPV...
    fx, fy, fz = spherical2carthesian(t, p)
    sx, sy, sz = project2normal(fx, fy, fz, nx, ny, nz)
    newray = propagate(Ix, Iy, Iz, sx, sy, sz)
    r = rho * cos(t)
    return true, newray, r
end

function bilambertian_path_shader(Ix::T, Iy::T, Iz::T, vx::T, vy::T, vz::T, nx::T, ny::T, nz::T ; rho::Array{T}, tau::Array{T}) where {T<:AbstractFloat}
    # note that vx, vy, and vz are not used but only exist so that input arguments are of the same pattern as in other shader functions...
    # note also that the whole spectrum undergoes one treatment: either reflect or transmit (obviously.)
    l = rho + tau
    r = 0.0
    if rand() > (sum(rho)/sum(l))
        nx, ny, nz = -1.0*nx, -1.0*ny, -1.0*nz
        r = tau
    else
        r = rho
    end
    t, p = sample_f() # could be replaced with a function RPVSampler if we knew how to implement importance sampling for RPV...
    fx, fy, fz = spherical2carthesian(t, p)
    sx, sy, sz = project2normal(fx, fy, fz, nx, ny, nz)
    newray = propagate(Ix, Iy, Iz, sx, sy, sz)
    r = r * cos(t)
    return true, newray, r
end

function rpv_path_shader(Ix::T, Iy::T, Iz::T, vx::T, vy::T, vz::T, nx::T, ny::T, nz::T ; rho_zero::Array{T}, rho_c::Array{T}, captheta::Array{T}, k::Array{T}) where {T<:AbstractFloat}
    t, p = sample_f() # could be replaced with a function RPVSampler if we knew how to implement importance sampling for RPV...
    fx, fy, fz = spherical2carthesian(t, p)
    sx, sy, sz = project2normal(fx, fy, fz, nx, ny, nz)
    phi = relazi(vx, vy, vz, nx, ny, nz, sx, sy, sz)
    theta = acos(vx*nx + vy*ny + vz*nz)
    theta_zero = acos(sx*nx + sy*ny + sz*nz)
    rho_sfc = rpv(theta_zero, theta, phi, rho_zero, rho_c, captheta, k)
    newray = propagate(Ix, Iy, Iz, sx, sy, sz)
    return true, newray, rho_sfc
end

function pure_reflection_path_shader(Ix::T, Iy::T, Iz::T, vx::T, vy::T, vz::T, nx::T, ny::T, nz::T ; rho::Array{T}) where {T<:AbstractFloat}
    sx, sy, sz = pure_reflection(vx, vy, vz, nx, ny, nz)
    newray = propagate(Ix, Iy, Iz, sx, sy, sz)
    return true, newray, rho
end

function lightsource_path_shader(Ix::T, Iy::T, Iz::T, vx::T, vy::T, vz::T, nx::T, ny::T, nz::T ; emission::Array{T}) where {T<:AbstractFloat}
    return false, Ray(Ix, Iy, Iz, vx, vy, vz), emission
end




