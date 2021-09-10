#########################################################################################################
#
#  The following shader functions have inputs that are all in the form (vx, vy, vz, nx, ny, nz, args...)
#  and return boolean, new ray, a multiplication/scattering factor (scalar or array)
#
#########################################################################################################



##############################################################################################
#
# Photon shaders... (single waveband, includes a probability to be absorbed)
#
##############################################################################################

function lambertian_full_fate_shader(Ix::T, Iy::T, Iz::T, vx::T, vy::T, vz::T, nx::T, ny::T, nz::T, rho::T) where {T<:AbstractFloat}
    # note that vx, vy, and vz are not used but only exist so that input arguments are of the same pattern as in other shader functions...
    scattering = false
    newray = Ray(Inf, Inf, Inf, Inf, Inf, Inf)
    r = 0.0
    if rand() <= rho
        scattering = true
        t, p = sample_f() # could be replaced with a function RPVSampler if we knew how to implement importance sampling for RPV...
	fx, fy, fz = spherical2carthesian(t, p)
	sx, sy, sz = project2normal(fx, fy, fz, nx, ny, nz)
	newray = propagate(Ix, Iy, Iz, sx, sy, sz)
	r = rho * cos(t)
    end
    return scattering, newray, r
end

function bilambertian_full_fate_shader(Ix::T, Iy::T, Iz::T, vx::T, vy::T, vz::T, nx::T, ny::T, nz::T, rho::T, tau::T) where {T<:AbstractFloat}
    # note that vx, vy, and vz are not used but only exist so that input arguments are of the same pattern as in other shader functions...
    scattering = false
    newray = Ray(Inf, Inf, Inf, Inf, Inf, Inf)
    r = 0.0
    omega = rho + tau # total scatter...
    if rand() <= omega
        scattering = true
        # choose either reflection or transmission...
        if rand() > (rho/omega)
            nx, ny, nz = -1.0*nx, -1.0*ny, -1.0*nz
        end
        # now, treats as BSDF and the total scatter is assigned to the direction of propagation
        t, p = sample_f()
        fx, fy, fz = spherical2carthesian(t, p)
        sx, sy, sz = project2normal(fx, fy, fz, nx, ny, nz)
        newray = propagate(Ix, Iy, Iz, sx, sy, sz)
        r = omega * cos(t)
    end
    return scattering, newray, r
end

function rpv_full_fate_shader(vx::T, vy::T, vz::T, nx::T, ny::T, nz::T, rho_zero::T, rho_c::T, captheta::T, k::T) where {T<:AbstractFloat}
    # note that vx, vy, and vz are not used but only exist so that input arguments are of the same pattern as in other shader functions...
    scattering = false
    newray = Ray(Inf, Inf, Inf, Inf, Inf, Inf)
    rho_sfc = 0.0
    if rand() <= rho_zero # 7SEPT2021: I PICKED RHO_ZERO BUT MAYBE IT SHOULD BE RHO_C INSTEAD?
        scattering = true
        t, p = sample_f() # could be replaced with a function RPVSampler if we knew how to implement importance sampling for RPV...
        fx, fy, fz = spherical2carthesian(t, p)
        vx, vy, vz = normalize(vx, vy, vz)
        nx, ny, nz = normalize(nx, ny, nz)
        sx, sy, sz = project2world(fx, fy, fz, nx, ny, nz)
        phi = relazi(vn, vy, vz, nx, ny, nz, sx, sy, sz)
        theta = abs(acos(vx*nx + vy*ny + vz*nz))
        theta_zero = abs(acos(sx*nx + sy*ny + sz*nz))
        rho_sfc = rpv(theta_zero, theta, phi, rho_zero, rho_c, captheta, k)
        newray = propagate(Ix, Iy, Iz, sx, sy, sz)
    end
    return scattering, newray, rho_sfc
end

function rpv_full_fate_shader(Ix::T, Iy::T, Iz::T, vx::T, vy::T, vz::T, nx::T, ny::T, nz::T, rho_zero::Array{T}, rho_c::Array{T}, captheta::Array{T}, k::Array{T}) where {T<:AbstractFloat}
    # note that vx, vy, and vz are not used but only exist so that input arguments are of the same pattern as in other shader functions...
    scattering = false
    newray = Ray(Inf, Inf, Inf, Inf, Inf, Inf)
    rho_sfc = similar(rho_zero)
    if rand() <= sum(rho_zero)/length(rho_zero) # 7SEPT2021: I PICKED RHO_ZERO BUT MAYBE IT SHOULD BE RHO_C INSTEAD?
        scattering = true
        t, p = sample_f() # could be replaced with a function RPVSampler if we knew how to implement importance sampling for RPV...
        fx, fy, fz = spherical2carthesian(t, p)
        vx, vy, vz = normalize(vx, vy, vz)
        nx, ny, nz = normalize(nx, ny, nz)
        sx, sy, sz = project2normal(fx, fy, fz, nx, ny, nz)
        phi = relazi(vx, vy, vz, nx, ny, nz, sx, sy, sz)
        theta = abs(acos(vx*nx + vy*ny + vz*nz))
        theta_zero = abs(acos(sx*nx + sy*ny + sz*nz))
        rho_sfc = rpv(theta_zero, theta, phi, rho_zero, rho_c, captheta, k)
        newray = propagate(Ix, Iy, Iz, sx, sy, sz)
    end
    return scattering, newray, rho_sfc
end


##############################################################################################
#
# Photon shaders... (spectra, includes probability of absorption)
#
##############################################################################################

function lambertian_full_fate_shader(Ix::T, Iy::T, Iz::T, vx::T, vy::T, vz::T, nx::T, ny::T, nz::T, rho::AbstractArray{T,1}) where {T<:AbstractFloat}
    # note that vx, vy, and vz are not used but only exist so that input arguments are of the same pattern as in other shader functions...
    scattering = false
    newray = Ray(Inf, Inf, Inf, Inf, Inf, Inf)
    r = 0.0
    if rand() <= sum(rho)/length(rho)
        scattering = true
        t, p = sample_f() # could be replaced with a function RPVSampler if we knew how to implement importance sampling for RPV...
	fx, fy, fz = spherical2carthesian(t, p)
	sx, sy, sz = project2normal(fx, fy, fz, nx, ny, nz)
	newray = propagate(Ix, Iy, Iz, sx, sy, sz)
	r = rho * cos(t)
    end
    return scattering, newray, r
end


function bilambertian_full_fate_shader(Ix::T, Iy::T, Iz::T, vx::T, vy::T, vz::T, nx::T, ny::T, nz::T, rho::AbstractArray{T,1}, tau::AbstractArray{T,1}) where {T<:AbstractFloat}
    # note that vx, vy, and vz are not used but only exist so that input arguments are of the same pattern as in other shader functions...
    scattering = false
    newray = Ray(Inf, Inf, Inf, Inf, Inf, Inf)
    r = 0.0
    omega = rho + tau # total scatter...
    if rand() <= sum(omega)/length(omega)
        scattering = true
        # choose either reflection or transmission...
        if rand() > (sum(rho)/sum(omega))
            nx, ny, nz = -1.0*nx, -1.0*ny, -1.0*nz
        end
        # now, treat as BSDF: the total scatter is assigned to the direction of propagation
        t, p = sample_f()
        fx, fy, fz = spherical2carthesian(t, p)
        sx, sy, sz = project2normal(fx, fy, fz, nx, ny, nz)
        newray = propagate(Ix, Iy, Iz, sx, sy, sz)
        r = omega * cos(t)
    end
    return scattering, newray, r
end



##############################################################################################
#
# Single-waveband shaders...
#
##############################################################################################

function lambertian_shader(Ix::T, Iy::T, Iz::T, vx::T, vy::T, vz::T, nx::T, ny::T, nz::T, rho::T) where {T<:AbstractFloat}
    # note that vx, vy, and vz are not used but only exist so that input arguments are of the same pattern as in other shader functions...
    t, p = sample_f() # could be replaced with a function RPVSampler if we knew how to implement importance sampling for RPV...
    fx, fy, fz = spherical2carthesian(t, p)
    sx, sy, sz = project2normal(fx, fy, fz, nx, ny, nz)
    newray = propagate(Ix, Iy, Iz, sx, sy, sz)
    r = rho * cos(t)
    return true, newray, r
end

function bilambertian_shader(Ix::T, Iy::T, Iz::T, vx::T, vy::T, vz::T, nx::T, ny::T, nz::T, rho::T, tau::T) where {T<:AbstractFloat}
    # note that vx, vy, and vz are not used but only exist so that input arguments are of the same pattern as in other shader functions...
    omega = rho + tau # total scatter...
    if rand() > (rho/omega)
        nx, ny, nz = -1.0*nx, -1.0*ny, -1.0*nz
    end
    # now, treats as BSDF and the total scatter is assigned to the direction of propagation
    t, p = sample_f()
    fx, fy, fz = spherical2carthesian(t, p)
    sx, sy, sz = project2normal(fx, fy, fz, nx, ny, nz)
    newray = propagate(Ix, Iy, Iz, sx, sy, sz)
    r = omega * cos(t)
    return true, newray, r
end

function pure_reflection_shader(vx::T, vy::T, vz::T, nx::T, ny::T, nz::T, rho::T) where {T<:AbstractFloat}
    sx, sy, sz = pure_reflection(vx, vy, vz, nx, ny, nz)
    newray = propagate(Ix, Iy, Iz, sx, sy, sz)
    return true, newray, rho
end

function rpv_shader(vx::T, vy::T, vz::T, nx::T, ny::T, nz::T, rho_zero::T, rho_c::T, captheta::T, k::T) where {T<:AbstractFloat}
    t, p = sample_f() # could be replaced with a function RPVSampler if we knew how to implement importance sampling for RPV...
    fx, fy, fz = spherical2carthesian(t, p)
    vx, vy, vz = normalize(vx, vy, vz)
    nx, ny, nz = normalize(nx, ny, nz)
    sx, sy, sz = project2world(fx, fy, fz, nx, ny, nz)
    phi = relazi(vn, vy, vz, nx, ny, nz, sx, sy, sz)
    theta = abs(acos(vx*nx + vy*ny + vz*nz))
    theta_zero = abs(acos(sx*nx + sy*ny + sz*nz))
    rho_sfc = rpv(theta_zero, theta, phi, rho_zero, rho_c, captheta, k)
    newray = propagate(Ix, Iy, Iz, sx, sy, sz)
    return true, newray, rho_sfc
end



###############################################################
#
#  Spectral shaders... (i.e. acting on spectra)
#
###############################################################

function lambertian_shader(Ix::T, Iy::T, Iz::T, vx::T, vy::T, vz::T, nx::T, ny::T, nz::T, rho::Array{T}) where {T<:AbstractFloat}
    # note that vx, vy, and vz are not used but only exist so that input arguments are of the same pattern as in other shader functions...
    t, p = sample_f() # could be replaced with a function RPVSampler if we knew how to implement importance sampling for RPV...
    fx, fy, fz = spherical2carthesian(t, p)
    sx, sy, sz = project2normal(fx, fy, fz, nx, ny, nz)
    newray = propagate(Ix, Iy, Iz, sx, sy, sz)
    r = rho * cos(t)
    return true, newray, r
end

function bilambertian_shader(Ix::T, Iy::T, Iz::T, vx::T, vy::T, vz::T, nx::T, ny::T, nz::T, rho::Array{T}, tau::Array{T}) where {T<:AbstractFloat}
    # note that vx, vy, and vz are not used but only exist so that input arguments are of the same pattern as in other shader functions...
    # note also that the whole spectrum undergoes one treatment: either reflect or transmit (obviously.)
    omega = rho + tau
    if rand() > (sum(rho)/sum(omega))
        nx, ny, nz = -1.0*nx, -1.0*ny, -1.0*nz
    end
    # treat as BSDF. All radiance goes to new direction...
    t, p = sample_f() # could be replaced with a function RPVSampler if we knew how to implement importance sampling for RPV...
    fx, fy, fz = spherical2carthesian(t, p)
    sx, sy, sz = project2normal(fx, fy, fz, nx, ny, nz)
    newray = propagate(Ix, Iy, Iz, sx, sy, sz)
    r = omega * cos(t)
    return true, newray, r
end

function rpv_shader(Ix::T, Iy::T, Iz::T, vx::T, vy::T, vz::T, nx::T, ny::T, nz::T, rho_zero::Array{T}, rho_c::Array{T}, captheta::Array{T}, k::Array{T}) where {T<:AbstractFloat}
    t, p = sample_f() # could be replaced with a function RPVSampler if we knew how to implement importance sampling for RPV...
    fx, fy, fz = spherical2carthesian(t, p)
    vx, vy, vz = normalize(vx, vy, vz)
    nx, ny, nz = normalize(nx, ny, nz)
    sx, sy, sz = project2normal(fx, fy, fz, nx, ny, nz)
    phi = relazi(vx, vy, vz, nx, ny, nz, sx, sy, sz)
    theta = abs(acos(vx*nx + vy*ny + vz*nz))
    theta_zero = abs(acos(sx*nx + sy*ny + sz*nz))
    rho_sfc = rpv(theta_zero, theta, phi, rho_zero, rho_c, captheta, k)
    newray = propagate(Ix, Iy, Iz, sx, sy, sz)
    return true, newray, rho_sfc
end

function pure_reflection_shader(Ix::T, Iy::T, Iz::T, vx::T, vy::T, vz::T, nx::T, ny::T, nz::T, rho::Array{T}) where {T<:AbstractFloat}
    sx, sy, sz = pure_reflection(vx, vy, vz, nx, ny, nz)
    newray = propagate(Ix, Iy, Iz, sx, sy, sz)
    return true, newray, rho
end

function lightsource_shader(Ix::T, Iy::T, Iz::T, vx::T, vy::T, vz::T, nx::T, ny::T, nz::T, emission::Array{T}) where {T<:AbstractFloat}
    return false, Ray(Ix, Iy, Iz, vx, vy, vz), emission
end


