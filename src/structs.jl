
import LinearAlgebra:normalize

function normalize(x::T, y::T, z::T) where {T<:Real}
  l::T = sqrt( *(x, x) + *(y, y) + *(z, z) )
  xn::T = x / l
  yn::T = y / l
  zn::T = z / l
  return xn, yn, zn
end


struct Asset{T<:AbstractFloat, S<:Integer}
    xoff::T
    yoff::T
    zoff::T
    beta::T
    gamma::T
    oidx::S
    pidx::S
end


struct Bvh{T<:AbstractFloat, S<:Integer, A<:AbstractArray, E<:AbstractArray}
    xmin::A
    ymin::A
    zmin::A
    xmax::A
    ymax::A
    zmax::A
    bvh2bb::E
    left_child::E
    right_child::E
    Bvh{T,S,A,E}(a::AbstractArray{T,1}, b::AbstractArray{T,1}, c::AbstractArray{T,1}, d::AbstractArray{T,1}, e::AbstractArray{T,1}, f::AbstractArray{T,1}, p::AbstractArray{S}, q::AbstractArray{S}, r::AbstractArray{S}) where {T,S,A,E} = new(a,b,c,d,e,f,p,q,r)
end
Bvh(a::AbstractArray, b::AbstractArray, c::AbstractArray, d::AbstractArray, e::AbstractArray, f::AbstractArray, p::AbstractArray, q::AbstractArray, r::AbstractArray) = Bvh{eltype(a), eltype(p), typeof(a), typeof(p)}(a,b,c,d,e,f,p,q,r)


struct Camera{T<:AbstractFloat, S<:Integer, A<:AbstractArray{T,1}}
  eyePoint::A
  lookAtPoint::A
  fov::T
  xResolution::S
  yResolution::S
  nBands::S
  nRaysPerPixel::S
  maxNumberOfCycles::S

  # derived attributes
  viewDirection::A
  u::A
  w::A
  aspectRatio::T
  viewPlaneHalfHeight::T
  viewPlaneHalfWidth::T
  viewPlaneBottomLeftPoint::A
  xIncVector::A
  yIncVector::A

  function Camera{T,S,A}(eyePoint::A, lookAtPoint::A, fov::T, xResolution::S, yResolution::S, nBands::S, nRaysPerPixel::S, maxNumberOfCycles::S) where {T,S,A}
    viewDirection = LinearAlgebra.normalize(lookAtPoint .- eyePoint)
    u = LinearAlgebra.normalize(cross(LinearAlgebra.normalize([0.000000, 0.000001, 1.0]), viewDirection))
    w = LinearAlgebra.normalize(cross(u, viewDirection))
    viewPlaneHalfWidth = tan(fov/2.0)
    aspectRatio = float(yResolution) / float(xResolution)
    viewPlaneHalfHeight = aspectRatio .* viewPlaneHalfWidth
    viewPlaneBottomLeftPoint = (eyePoint .+ viewDirection) .- (w .* viewPlaneHalfHeight) .- (u .* viewPlaneHalfWidth)
    xIncVector = (u .* 2.0 .* viewPlaneHalfWidth) ./ float(xResolution)
    yIncVector = (w .* 2.0 .* viewPlaneHalfHeight) ./ float(yResolution)
    new(eyePoint, lookAtPoint, fov, xResolution, yResolution, nBands, nRaysPerPixel, maxNumberOfCycles, viewDirection, u, w, aspectRatio, viewPlaneHalfHeight, viewPlaneHalfWidth, viewPlaneBottomLeftPoint, xIncVector, yIncVector)
  end
end
Camera(v::A, w::A, p::T, q::S, r::S, s::S, t::S, u::S) where {T<:AbstractFloat, S<:Integer, A<:AbstractArray} = Camera{T,S,A}(v,w,p,q,r,s,t,u)


struct Coord{S<:Integer}
    x::S
    y::S
end


abstract type Geometry{T<:AbstractFloat, S<:Integer, A<:AbstractArray, E<:AbstractArray} end


struct Ball{T,S,A,E} <: Geometry{T,S,A,E}
    p1x::A
    p1y::A
    p1z::A
    r::A
    rsq::A
    mtl::E
    function Ball{T,S,A,E}(p::AbstractArray{T,2}, radius::AbstractArray{T,1}, mtl::AbstractArray{S,1}) where {T,S,A,E}
        n_balls = size(p, 1)
        
        p1x::Array{T, 1} = p[:,1]
        p1y::Array{T, 1} = p[:,2]
        p1z::Array{T, 1} = p[:,3]
        
        r::Array{T, 1} = radius
        rsq::Array{T, 1} = radius.^2
        
        mtl::Array{S, 1} = mtl

        new(p1x, p1y, p1z, r, rsq, mtl)
    end
end
Ball(p::AbstractArray{T,2}, radius::AbstractArray{T,1}, mtl::AbstractArray{S,1}) where {T<:AbstractFloat, S<:Integer} = Ball{eltype(p), eltype(mtl), Array{T,1}, Array{S,1}}(p,radius,mtl)


struct Disk{T,S,A,E} <: Geometry{T,S,A,E}
    p1x::A
    p1y::A
    p1z::A
    nx::A
    ny::A
    nz::A
    r::A
    rsq::A   
    
    mtl::E
    
    function Disk{T,S,A,E}(p::AbstractArray{T,2}, n::AbstractArray{T,2}, radius::AbstractArray{T,1}, mtl::AbstractArray{S,1}) where {T,S,A,E}
        n_disks = size(p, 1)
        
        p1x::Array{T, 1} = p[:,1]
        p1y::Array{T, 1} = p[:,2]
        p1z::Array{T, 1} = p[:,3]
        
        nx::Array{T, 1} = n[:,1]
        ny::Array{T, 1} = n[:,2]
        nz::Array{T, 1} = n[:,3]
        
        r::Array{T, 1} = radius
        rsq::Array{T, 1} = radius.^2
        
        mtl::Array{S, 1} = mtl

        new(p1x, p1y, p1z, nx, ny, nz, r, rsq, mtl)
    end
end
Disk(p::AbstractArray{T,2}, n::AbstractArray{T,2}, radius::AbstractArray{T,1}, mtl::AbstractArray{S,1}) where {T<:AbstractFloat, S<:Integer} = Disk{eltype(p), eltype(mtl), Array{T,1}, Array{S,1}}(p,n,radius,mtl)


struct Cylinder{T,S,A,E} <: Geometry{T,S,A,E}
    m1x::A
    m1y::A
    m1z::A
    m2x::A
    m2y::A
    m2z::A
    radius::A
    u::A
    v::A
    w::A    
    length::A
    beta::A
    gamma::A
    rsq::A
    
    mtl::E
    
    function Cylinder{T,S,A,E}(m::AbstractArray{T,2}, l::AbstractArray{S,2}, radius::AbstractArray{T,1}, mtl::AbstractArray{S,1}) where {T,S,A,E}
        n_cylinders = size(l, 1)
        
        m1x::Array{T, 1} = zeros(n_cylinders)
        m1y::Array{T, 1} = zeros(n_cylinders)
        m1z::Array{T, 1} = zeros(n_cylinders)
        
        m2x::Array{T, 1} = zeros(n_cylinders)
        m2y::Array{T, 1} = zeros(n_cylinders)
        m2z::Array{T, 1} = zeros(n_cylinders)
        
        radius::Array{T, 1} = radius        
        u::Array{T, 1} = zeros(n_cylinders)
        v::Array{T, 1} = zeros(n_cylinders)
        w::Array{T, 1} = zeros(n_cylinders)
        length::Array{T, 1} = zeros(n_cylinders)
        beta::Array{T, 1} = zeros(n_cylinders)
        gamma::Array{T, 1} = zeros(n_cylinders)
        rsq::Array{T, 1} = zeros(n_cylinders)
        
        mtl::Array{S, 1} = mtl
        
        for i = 1:n_cylinders
            idx1 = l[i,1]; idx2 = l[i,2]
            m1x[i] = m[idx1,1]; m1y[i] = m[idx1,2]; m1z[i] = m[idx1,3]
            m2x[i] = m[idx2,1]; m2y[i] = m[idx2,2]; m2z[i] = m[idx2,3]
            dx = m[idx2,1] - m[idx1,1]
            dy = m[idx2,2] - m[idx1,2]
            dz = m[idx2,3] - m[idx1,3]
            length_ = sqrt(dx^2 + dy^2 + dz^2)
            length[i] = length_
            beta[i] = acos(dz / length_)
            gamma[i] = atan(dy, dx) # formerly atan2
            rsq[i] = radius[i]^2
            # compute a 'normal' vector that is pointing up
            q = [dx, dy, dz]
            nv = normalize(cross([0.0034, 0.0071, 1.0], q))
            vq = cross(nv, q)
            u[i] = vq[1]
            v[i] = vq[2]
            w[i] = vq[3]            
        end
        new(m1x, m1y, m1z, m2x, m2y, m2z, radius, u, v, w, length, beta, gamma, rsq, mtl)
    end
end
Cylinder(m::AbstractArray{T,2}, l::AbstractArray{S,2}, radius::AbstractArray{T,1}, mtl::AbstractArray{S,1}) where {T<:AbstractFloat, S<:Integer} = Cylinder{eltype(m), eltype(l), Array{T,1}, Array{S,1}}(m,l,radius,mtl)


struct Cone{T,S,A,E} <: Geometry{T,S,A,E}
    m1x::A
    m1y::A
    m1z::A    
    m2x::A
    m2y::A
    m2z::A 
    radius_base::A
    radius_tip::A
    u::A
    v::A
    w::A
    beta::A
    gamma::A
    length::A
    a::A
    aa::A
    r1a::A
    r1sq::A
    
    mtl::E
    
    function Cone{T,S,A,E}(m::AbstractArray{T,2}, l::AbstractArray{S,2}, radius_base::AbstractArray{T,1}, radius_tip::AbstractArray{T,1}, mtl::AbstractArray{S,1}) where {T,S,A,E}
        n_cones = size(l, 1)
        
        m1x::Array{T, 1} = zeros(n_cones)
        m1y::Array{T, 1} = zeros(n_cones)
        m1z::Array{T, 1} = zeros(n_cones)
        
        m2x::Array{T, 1} = zeros(n_cones)
        m2y::Array{T, 1} = zeros(n_cones)
        m2z::Array{T, 1} = zeros(n_cones)
        
        radius_base::Array{T, 1} = radius_base
        radius_tip::Array{T, 1} = radius_tip
        
        u::Array{T, 1} = zeros(n_cones)
        v::Array{T, 1} = zeros(n_cones)
        w::Array{T, 1} = zeros(n_cones)
        beta::Array{T, 1} = zeros(n_cones)
        gamma::Array{T, 1} = zeros(n_cones)
        length::Array{T, 1} = zeros(n_cones)
        a::Array{T, 1} = zeros(n_cones)
        aa::Array{T, 1} = zeros(n_cones)
        r1a::Array{T, 1} = zeros(n_cones)
        r1sq::Array{T, 1} = zeros(n_cones)
        
        mtl::Array{S, 1} = mtl
        
        for i = 1:n_cones
            idx1 = l[i,1]; idx2 = l[i,2]            
            m1x[i] = m[idx1,1]; m1y[i] = m[idx1,2]; m1z[i] = m[idx1,3]
            m2x[i] = m[idx2,1]; m2y[i] = m[idx2,2]; m2z[i] = m[idx2,3]
            dx = m[idx2,1] - m[idx1,1]
            dy = m[idx2,2] - m[idx1,2]
            dz = m[idx2,3] - m[idx1,3]
            length_ = sqrt(dx^2 + dy^2 + dz^2)
            length[i] = length_
            beta[i] = acos(dz / length_)
            gamma[i] = atan(dy, dx) # formerly atan2
            a_ = (radius_tip[i] - radius_base[i]) / length_ # taper of the stem
            a[i] = a_
            aa[i] = a_^2
            r1a[i] = radius_base[i]*a_
            r1sq[i] = radius_base[i]^2
            # compute a 'normal' vector that is pointing up
            q = [dx, dy, dz]
            nv = normalize(cross([0.0034, 0.0071, 1.0], q))
            vq = cross(nv, q)
            u[i] = vq[1]
            v[i] = vq[2]
            w[i] = vq[3]            
        end
        new(m1x, m1y, m1z, m2x, m2y, m2z, radius_base, radius_tip, u, v, w, beta, gamma, length, a, aa, r1a, r1sq, mtl)
    end
end
Cone(m::AbstractArray{T,2}, l::AbstractArray{S,2}, radius_base::AbstractArray{T,1}, radius_tip::AbstractArray{T,1}, mtl::AbstractArray{S,1}) where {T<:AbstractFloat, S<:Integer} = Cone{eltype(m), eltype(l), Array{T,1}, Array{S,1}}(m,l,radius_base,radius_tip,mtl)


struct Mesh{T,S,A,E} <: Geometry{T,S,A,E}
    v1x::A
    v1y::A
    v1z::A

    v2x::A
    v2y::A
    v2z::A

    v3x::A
    v3y::A
    v3z::A

    mtl::E

    ux::A
    uy::A
    uz::A

    wx::A
    wy::A
    wz::A

    nx::A
    ny::A
    nz::A

    uu::A
    uw::A
    ww::A
    D::A

    function Mesh{T,S,A,E}(v::AbstractArray{T, 2}, t::AbstractArray{S, 2}, m::AbstractArray{S, 1}) where {T,S,A,E}
        n_triangles = size(t, 1)
        v1x::Array{T, 1} = zeros(n_triangles)
        v1y::Array{T, 1} = zeros(n_triangles)
        v1z::Array{T, 1} = zeros(n_triangles)
        v2x::Array{T, 1} = zeros(n_triangles)
        v2y::Array{T, 1} = zeros(n_triangles)
        v2z::Array{T, 1} = zeros(n_triangles)
        v3x::Array{T, 1} = zeros(n_triangles)
        v3y::Array{T, 1} = zeros(n_triangles)
        v3z::Array{T, 1} = zeros(n_triangles)

        mtl::Array{S, 1} = m

        ux::Array{T, 1} = zeros(n_triangles)
        uy::Array{T, 1} = zeros(n_triangles)
        uz::Array{T, 1} = zeros(n_triangles)

        wx::Array{T, 1} = zeros(n_triangles)
        wy::Array{T, 1} = zeros(n_triangles)
        wz::Array{T, 1} = zeros(n_triangles)

        nx::Array{T, 1} = zeros(n_triangles)
        ny::Array{T, 1} = zeros(n_triangles)
        nz::Array{T, 1} = zeros(n_triangles)

        uu::Array{T, 1} = zeros(n_triangles)
        uw::Array{T, 1} = zeros(n_triangles)
        ww::Array{T, 1} = zeros(n_triangles)
        D::Array{T, 1} = zeros(n_triangles)
        for i = 1:n_triangles
            idx1 = t[i,1]; idx2 = t[i,2]; idx3 = t[i,3]
            v1x[i] = v[idx1,1]; v1y[i] = v[idx1,2]; v1z[i] = v[idx1,3]
            v2x[i] = v[idx2,1]; v2y[i] = v[idx2,2]; v2z[i] = v[idx2,3]
            v3x[i] = v[idx3,1]; v3y[i] = v[idx3,2]; v3z[i] = v[idx3,3]

            # precompute some variables...
            ux_ = v2x[i] - v1x[i]
            uy_ = v2y[i] - v1y[i]
            uz_ = v2z[i] - v1z[i]
            wx_ = v3x[i] - v1x[i]
            wy_ = v3y[i] - v1y[i]
            wz_ = v3z[i] - v1z[i]
            u_ = [ux_, uy_, uz_]
            w_ = [wx_, wy_, wz_]
            n_ = cross(u_, w_)
            n2_ = n_.^2
            sm_ = sum(n2_)
            l_ = sqrt(sm_)
            n_ = n_./l_
            uu_ = dot(u_, u_)
            uw_ = dot(u_, w_)
            ww_ = dot(w_, w_)
            uw2_ = (uw_ * uw_)
            uuww_ = (uu_ * ww_)
            D_ = uw2_ - uuww_

            ux[i] = ux_
            uy[i] = uy_
            uz[i] = uz_

            wx[i] = wx_
            wy[i] = wy_
            wz[i] = wz_

            nx[i] = n_[1]
            ny[i] = n_[2]
            nz[i] = n_[3]

            uu[i] = uu_
            uw[i] = uw_
            ww[i] = ww_
            D[i] = D_
        end
        new(v1x, v1y, v1z, v2x, v2y, v2z, v3x, v3y, v3z, mtl, ux, uy, uz, wx, wy, wz, nx, ny, nz, uu, uw, ww, D)
    end
end
Mesh(v::AbstractArray{T,2}, t::AbstractArray{S,2}, m::AbstractArray{S,1}) where {T<:AbstractFloat, S<:Integer} = Mesh{eltype(v), eltype(t), Array{T,1}, Array{S,1}}(v,t,m)


struct Point{T<:AbstractFloat}
    x::T
    y::T
    z::T
end


struct Ray{T<:AbstractFloat}
    x::T
    y::T
    z::T
    dx::T
    dy::T
    dz::T
end


struct Hook{T<:AbstractFloat, S<:Integer}
    xoff::T
    yoff::T
    zoff::T
    beta::T
    gamma::T
    idx::S
    Hook{T,S}(x::T, y::T, z::T, beta::T, gamma::T, idx::S) where {T,S} = new(x,y,z,beta,gamma,idx)
end
Hook(x::AbstractFloat, y::AbstractFloat, z::AbstractFloat, beta::AbstractFloat, gamma::AbstractFloat, idx::Integer) = Hook{typeof(x), typeof(idx)}(x,y,z,beta,gamma,idx)


