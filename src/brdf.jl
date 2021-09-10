function cross_scalars(ax::T, ay::T, az::T, bx::T, by::T, bz::T) where {T<:AbstractFloat}
    cx = ay*bz-az*by
    cy = az*bx-ax*bz
    cz = ax*by-ay*bx
    return cx, cy, cz
end

function normalize_scalars(x::T, y::T, z::T) where {T<:AbstractFloat}
    l = sqrt(x*x + y*y + z*z)
    x = x / l
    y = y / l
    z = z / l
    return x, y, z
end

function project2normal(fx::T, fy::T, fz::T, nx::T, ny::T, nz::T) where {T<:AbstractFloat}
    ax::T = 0.0034
    ay::T = 0.0071
    az::T = 1.0000
    bx, by, bz = cross_scalars(ax, ay, az, nx, ny, nz)
    cx, cy, cz = normalize_scalars(bx, by, bz)
    dx, dy, dz = cross_scalars(cx, cy, cz, nx, ny, nz)
    ex, ey, ez = normalize_scalars(dx, dy, dz)
    gx, gy, gz = normalize_scalars(fx, fy, fz)
    hx = gx * ex + gy * cx + gz * nx
    hy = gx * ey + gy * cy + gz * ny
    hz = gx * ez + gy * cz + gz * nz
    return hx, hy, hz
end

function sample_f()
    x1 = rand()
    x2 = rand()
    phi = x1*2.0*pi
    theta = acos(x2)
    return theta, phi
end

function propagate(Ix::T, Iy::T, Iz::T, dx::T, dy::T, dz::T) where {T<:AbstractFloat}
    x = Ix + dx * SMALL_NUM
    y = Iy + dy * SMALL_NUM
    z = Iz + dz * SMALL_NUM
    ray = Ray(x, y, z, dx, dy, dz)
    return ray
end

function compute_propagation_of_reflectance(Ix::T, Iy::T, Iz::T, nx::T, ny::T, nz::T) where {T<:AbstractFloat}
    t, p = sample_f()
    fx, fy, fz = spherical2carthesian(t, p)
    dx, dy, dz = project2normal(fx, fy, fz, nx, ny, nz)
    x = Ix + dx * SMALL_NUM
    y = Iy + dy * SMALL_NUM
    z = Iz + dz * SMALL_NUM
    ray = Ray(x, y, z, dx, dy, dz)
    return ray
end

function compute_propagation_of_transmittance(Ix::T, Iy::T, Iz::T, nx::T, ny::T, nz::T) where {T<:AbstractFloat}
    t, p = sample_f()
    fx, fy, fz = spherical2carthesian(t, p)
    dx, dy, dz = project2normal(fx, fy, fz, -nx, -ny, -nz)
    x = Ix + dx * SMALL_NUM
    y = Iy + dy * SMALL_NUM
    z = Iz + dz * SMALL_NUM
    ray = Ray(x, y, z, dx, dy, dz)
    return ray
end

function distr(nx::T, ny::T, nz::T, hx::T, hy::T, hz::T, alpha::T) where {T<:AbstractFloat}
    alpha2 = alpha^2
    NoH = nx*hx+ny*hy+nz*hz
    NoH2 = NoH*NoH
    den = NoH2 * alpha2 + (1 - NoH2)
    chi = NoH > 0 ? 1 : 0
    return (chi * alpha2) / (np.pi * den * den)
end

function geom(vx::T, vy::T, vz::T, nx::T, ny::T, nz::T, hx::T, hy::T, hz::T, alpha::T) where {T<:AbstractFloat}
    VoH = vx*hx + vy*hy + vz*hz
    VoN = vx*nx + vy*ny + vz*nz
    m = VoH/VoN
    chi = m > 0 ? 1 : 0
    VoH2 = VoH * VoH
    tan2 = (1 - VoH2) / VoH2
    return (chi*2) / (1 + sqrt(1 + alpha * alpha * tan2))
end

function fresnel(vx::T, vy::T, vz::T, hx::T, hy::T, hz::T, eta::T, k::T) where {T<:AbstractFloat}
    VoH = vx*hx + vy*hy + vz*hz
    k2 = k^2
    num = (eta-1)^2 + 4*eta*(1-VoH)^5 + k2
    den = (eta+1)^2 + k^2
    return num/den
end

function cooktorrance(vx::T, vy::T, vz::T, nx::T, ny::T, nz::T, sx::T, sy::T, sz::T, alpha::T, eta::T, k::T) where {T<:AbstractFloat}
    h = compute_half_vector(vx, vy, vz, sx, sy, sz)
    F = fresnel(vx, vh, vz, hx, hy, hz, eta, k)
    G = geom(vx, vy, vz, nx, ny, nz, hx, hy, hz, alpha) * geom(sx, sy, sz, nx, ny, nz, hx, hy, hz, alpha)
    D = distr(nx, ny, nz, hx, hy, hz, alpha)
    cosT = sx*nx + sy*ny + sz*nz
    sinT = sqrt(1-cosT*cosT)
    nom = F*G*D*sinT
    NoV = nx*vx + ny*vy + nz*vz
    HoN = nx*hx + ny*hy + nz*hz
    den = 4*NoV*HoN+0.05
    fs = nom/den
    return fs
end

function compute_half_vector(vx::T, vy::T, vz::T, sx::T, sy::T, sz::T) where {T<:AbstractFloat}
    hx, hy, hz = normalize_scalars(vx+sx, vy+sy, vz+sz)
    return hx, hy, hz
end

function ct_sampler(alpa::T) where {T<:AbstractFloat}
    jota1, jota2 = rand(2)
    num = alpha * sqrt(jota1)
    den = sqrt(1-jota1)
    theta = atan(num/den)
    phi = jota2*2.0*pi
    return theta, phi
end

function pure_reflection(vx::T, vy::T, vz::T, nx::T, ny::T, nz::T) where {T<:AbstractFloat}
    NoV = nx*vx + ny*vy + nz*vz
    rx = 2.0 * NoV * nx - vx
    ry = 2.0 * NoV * ny - vy
    rz = 2.0 * NoV * nz - vz
    return rx, ry, rz
end

function rpv(theta_zero::T, theta::T, phi::T, rho_zero::T, rho_c::T, captheta::T, k::T) where {T<:AbstractFloat}
    cos_theta = cos(theta)
    cos_theta_zero = cos(theta_zero)
    cos_phi = cos(phi)
    tan_theta_zero = tan(theta_zero)
    tan_theta = tan(theta)
    cosg = cos_theta*cos_theta_zero + sin(theta)*sin(theta_zero)*cos_phi
    captheta_square = captheta^2
    G = ((tan_theta_zero^2 + tan_theta^2) - (2*tan_theta_zero*tan_theta*cos_phi))^0.5
    M = (cos_theta_zero^(k-1) * cos_theta^(k-1)) / (cos_theta_zero+cos_theta)^(k-1)
    FHG = (1-captheta_square) / (1+2*captheta*cosg+captheta_square)^1.5
    H = 1+(1-rho_c)/(1+G)
    rho_sfc = rho_zero * M * FHG * H
    return rho_sfc
end

function rpv(theta_zero::T, theta::T, phi::T, rho_zero::Array{T}, rho_c::Array{T}, captheta::Array{T}, k::Array{T}) where {T<:AbstractFloat}
    cos_theta = cos(theta)
    cos_theta_zero = cos(theta_zero)
    cos_phi = cos(phi)
    tan_theta_zero = tan(theta_zero)
    tan_theta = tan(theta)
    cosg = cos_theta*cos_theta_zero + sin(theta)*sin(theta_zero)*cos_phi
    G = ((tan_theta_zero^2 + tan_theta^2) - (2*tan_theta_zero*tan_theta*cos_phi))^0.5
    n = size(k,1)
    rho_sfc = Array{T}(undef, n)
    for i=1:n
        captheta_square = captheta[i]^2
        M = (cos_theta_zero^(k[i]-1) * cos_theta^(k[i]-1)) / (cos_theta_zero+cos_theta)^(k[i]-1)
        FHG = (1 - captheta_square) / (1 + 2 * captheta[i] * cosg + captheta_square)^1.5
        H = 1 + (1 - rho_c[i]) / (1 + G)
        rho_sfc[i] = rho_zero[i] * M * FHG * H
    end
    return rho_sfc
end

function rpv(vx::T, vy::T, vz::T, nx::T, ny::T, nz::T, sx::T, sy::T, sz::T, rho_zero::T, rho_c::T, captheta::T, k::T) where {T<:AbstractFloat}
    phi = relazi(vn, vy, vz, nx, ny, nz, sx, sy, sz)
    theta = acos(vx*nx + vy*ny + vz*nz)
    theta_zero = acos(sx*nx + sy*ny + sz*nz)
    rho_sfc = RPV(theta_zero, theta, phi, rho_zero, rho_c, captheta, k)
    return rho_sfc
end


