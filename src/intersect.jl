function hit_box(minx::T, miny::T, minz::T, maxx::T, maxy::T, maxz::T, ray::Ray{T}) where {T<:AbstractFloat}
    minx2v1 = minx - ray.x
    maxx2v1 = maxx - ray.x
    miny2v1 = miny - ray.y
    maxy2v1 = maxy - ray.y
    minz2v1 = minz - ray.z
    maxz2v1 = maxz - ray.z

    tmin = minx2v1 / ray.dx
    tmax = maxx2v1 / ray.dx
    if tmin > tmax
        tmin, tmax = tmax, tmin
    end

    tymin = miny2v1 / ray.dy
    tymax = maxy2v1 / ray.dy
    if tymin > tymax
        tymin, tymax = tymax, tymin
    end

    if (tmin > tymax) || (tymin > tmax)
        return false
    end

    if tymin > tmin
        tmin = tymin
    end

    if tymax < tmax
        tmax = tymax
    end

    tzmin = minz2v1 / ray.dz
    tzmax = maxz2v1 / ray.dz
    if tzmin > tzmax
        tzmin, tzmax = tzmax, tzmin
    end

    if (tmin > tzmax) || (tzmin > tmax)
        return false
    end

    return true
end

function intersect_walls(walls::Walls, ray::Ray)
    r1::Float64 = Inf
    r2::Float64 = Inf
    r3::Float64 = Inf

    xwall::Char = 'o'
    ywall::Char = 'o'
    zwall::Char = 'o'
    wall::Char = 'x'

    if ray.dx > 0.0
        r1 = (walls.east - ray.x) / ray.dx
        xwall = 'e'
    elseif ray.dx < 0.0
        r1 = (walls.west - ray.x) / ray.dx
        xwall = 'w'
    end

    if ray.dy > 0.0
        r2 = (walls.north - ray.y) / ray.dy
        ywall = 'n'
    elseif ray.dy < 0.0
        r2 = (walls.south - ray.y) / ray.dy
        ywall = 's'
    end

    if ray.dz > 0.0
        r3 = (walls.top - ray.z) / ray.dz
        zwall = 't'
    elseif ray.dz < 0.0
        r3 = (walls.bottom - ray.z) / ray.dz
        zwall = 'b'
    end

    r = min(r1, r2, r3)

    if r === r1
        wall = xwall
    elseif r === r2
        wall = ywall
    else
        wall = zwall
    end

    x = ray.x + ray.dx * r
    y = ray.y + ray.dy * r
    z = ray.z + ray.dz * r

    I::Point = Point(x, y, z)

    return I, wall, r
end


function cycle_ray(walls::Walls, ray::Ray)
    # How to type multiple return values?
    I, wall, r = intersect_walls(walls, ray)

    x::Float64 = I.x
    y::Float64 = I.y
    z::Float64 = I.z

    xsame::Bool = false
    ysame::Bool = false
    zsame::Bool = false

    # make these constants?
    small_num::Float64 = 0.000001
    u::Float64 = 100.0

    if I.x <= (walls.west + small_num)
        x = I.x + walls.xsize
    elseif I.x >= (walls.east - small_num)
        x = I.x - walls.xsize
    else
        xsame = true
    end

    if I.y <= (walls.south + small_num)
        y = I.y + walls.ysize
    elseif I.y >= (walls.north - small_num)
        y = I.y - walls.ysize
    else
        ysame = true
    end

    if I.z <= (walls.bottom + small_num)
        z = I.z + walls.zsize
    elseif I.z >= (walls.top - small_num)
        z = I.z - walls.zsize
    else
        zsame = true
    end

    # assert that at least one booleans evaluates to true
    # skipped for now...
    new_ray::Ray = Ray(x, y, z, ray.dx, ray.dy, ray.dz)

    return new_ray
end


function intersect_bvh(bvh::Bvh{T,S,A,E}, ray::Ray{T}) where {T<:AbstractFloat, S<:Integer, A<:AbstractArray{T,1}, E<:AbstractArray{S,1}}
    #NOTE:these two vector take about a third of the time of this function and furthermore the pushes take a long time too...
    # the hit_box line takes surprisingly little, with only 586/4505
    q::Vector = [1]
    leaf_idxs::Vector{S} = []
    while size(q, 1) > 0
        i = pop!(q)
        minx, miny, minz, maxx, maxy, maxz = bvh.xmin[i], bvh.ymin[i], bvh.zmin[i], bvh.xmax[i], bvh.ymax[i], bvh.zmax[i]
        if hit_box(minx, miny, minz, maxx, maxy, maxz, ray)::Bool
            if bvh.left_child[i] > 1 && bvh.right_child[i] > 1
                q = push!(q, bvh.left_child[i], bvh.right_child[i])
            end
            if bvh.bvh2bb[i] > 0
                leaf_idxs = push!(leaf_idxs, bvh.bvh2bb[i])
            end
        end
    end
    return leaf_idxs
end


function solve_quadratic(a::T, b::T, c::T) where {T<:AbstractFloat}
    x0::T = 0.0
    x1::T = 0.0
    discr::T = b*b - 4*a*c
    if discr < 0
        return false, x0, x1
    elseif (discr==0)
        x0 = x1 = -0.5*b/a
    else
        q = b>0 ? -0.5*(b+sqrt(discr)) : -0.5*(b-sqrt(discr))
        x0 = q / a
        x1 = c / q
    end
    if x0 > x1
        x0, x1 = x1, x0
    end
    return true, x0, x1
end


function intersect(balls::Ball{T,S,A,E}, i::S, ray::Ray{T}) where {T<:AbstractFloat, S<:Integer, A<:AbstractArray, E<:AbstractArray}
    Lx::T = ray.x - balls.p1x[i]
    Ly::T = ray.y - balls.p1y[i]
    Lz::T = ray.z - balls.p1z[i]
    a::T = ray.dx*ray.dx + ray.dy*ray.dy + ray.dz*ray.dz
    b::T = 2 * (Lx*ray.dx + Ly*ray.dy + Lz*ray.dz)
    c::T = Lx*Lx + Ly*Ly + Lz*Lz - balls.rsq[i]
    hit, t0, t1 = solve_quadratic(a, b, c)
    if !hit
        nullpoint = Point{T}(Inf, Inf, Inf)
        return false, nullpoint, nullpoint, Inf
    end
    if t0 < 0 # if t0 is negative then use t1
	t0 = t1
    end
    if t0 >= 0
	Ix = ray.x + ray.dx * t0
	Iy = ray.y + ray.dy * t0
	Iz = ray.z + ray.dz * t0
	I = Point{T}(Ix, Iy, Iz)
	dx = Ix - balls.p1x[i]
        dy = Iy - balls.p1y[i]
        dz = Iz - balls.p1z[i]
        l = sqrt(dx^2 + dy^2 + dz^2)
        dx /= l
        dy /= l
        dz /= l
        normal::Point = Point{T}(dx, dy, dz)
	return hit, I, normal, t0
    end
    nullpoint = Point{T}(Inf, Inf, Inf)
    return false, nullpoint, nullpoint, Inf
end


function intersect(disks::Disk{T,S,A,E}, i::S, ray::Ray{T}) where {T<:AbstractFloat, S<:Integer, A<:AbstractArray, E<:AbstractArray}
    denom::T = disks.nx[i]*ray.dx + disks.ny[i]*ray.dy + disks.nz[i]*ray.dz
    if abs(denom) > SMALL_NUM
        # check point of intersection on the plane...
        p0l0x = disks.p1x[i] - ray.x
        p0l0y = disks.p1y[i] - ray.y
        p0l0z = disks.p1z[i] - ray.z
        t = (p0l0x*disks.nx[i] + p0l0y*disks.ny[i] + p0l0z*disks.nz[i]) / denom
        if t < 0
            nullpoint = Point{T}(Inf, Inf, Inf)
            return false, nullpoint, nullpoint, Inf
        end
        Ix = ray.x + ray.dx * t
        Iy = ray.y + ray.dy * t
        Iz = ray.z + ray.dz * t
        # see if it is on the disk...
        vx = Ix - disks.p1x[i]
        vy = Iy - disks.p1y[i]
        vz = Iz - disks.p1z[i]
        vv = vx*vx + vy*vy + vz*vz
        hit = vv <= disks.rsq[i]
        I = Point{T}(Ix, Iy, Iz)
	normal::Point{T} = Point{T}(disks.nx[i], disks.ny[i], disks.nz[i])
        return hit, I, normal, t
    end
    nullpoint = Point{T}(Inf, Inf, Inf)
    return false, nullpoint, nullpoint, Inf
end


function intersect(cone::Cone{T,S,A,E}, i::S, ray::Ray) where {T<:AbstractFloat, S<:Integer, A<:AbstractArray, E<:AbstractArray}
  # intersect ray with cone/cylinder...
  xoff = cone.m1x[i]
  yoff = cone.m1y[i]
  zoff = cone.m1z[i]
  beta = cone.beta[i]
  gamma = cone.gamma[i]
  
  tray = transform(ray, -xoff, -yoff, -zoff, -beta, -gamma)

  a::Float64 = tray.dx^2 + tray.dy^2 - tray.dz^2 * cone.aa[i]
  b::Float64 = 2 * tray.x * tray.dx + 2 * tray.y * tray.dy - 2 * cone.r1a[i] * tray.dz - 2 * tray.z * ray.dz * cone.aa[i]
  c::Float64 = tray.x^2 + tray.y^2 - cone.r1sq[i] - 2 * cone.r1a[i] * tray.z - cone.aa[i] * tray.z^2
  d::Float64 = b^2 - (4*a*c)

  t1::Float64 = Inf
  t2::Float64 = Inf
  t::Float64 = Inf

  if d > 0
    t1 = (-1.0*b - sqrt(d))/(2*a)
    t2 = (-1.0*b + sqrt(d))/(2*a)
  end

  # choose the smallest, non-negative solution for t
  if t1 <= 0
    if t2 > 0
      t = t2
    end
  else
    if t2 <= 0
      t = t1
    elseif t1 < t2
      t = t1
    else
      t = t2
    end
  end
  # t1 and t2 can both be negative, hence t can be negative!

  hit::Bool = false
  tI = Point{T}(tray.x + tray.dx * t, tray.y + tray.dy * t, tray.z + tray.dz * t)
  if (t > 0) && (tI.z > 0.0) && (tI.z < cone.length[i])
    hit = true
  end

  I = Point{T}(ray.x + ray.dx * t, ray.y + ray.dy * t, ray.z + ray.dz * t)
  xx::Float64 = tI.x^2
  yy::Float64 = tI.y^2
  rcirc::Float64 = sqrt(xx + yy)
  nx::Float64, ny::Float64, nz::Float64 = normalize(tI.x, tI.y, cone.a[i] * -1.0 * rcirc) 
  nx, ny, nz = transform(nx, ny, nz, xoff, yoff, zoff, beta, gamma)
  normal::Point = Point{T}(nx, ny, nz)

  return hit, I, normal, t
end


function intersect(cylinder::Cylinder{T,S,A,E}, i::S, ray::Ray) where {T<:AbstractFloat, S<:Integer, A<:AbstractArray, E<:AbstractArray}
  # intersect ray with cylinder...
  
  xoff = cylinder.m1x[i]
  yoff = cylinder.m1y[i]
  zoff = cylinder.m1z[i]
  beta = cylinder.beta[i]
  gamma = cylinder.gamma[i]
    
  tray = transform(ray, -xoff, -yoff, -zoff, -beta, -gamma)
    
  a::Float64 = tray.dx^2 + tray.dy^2
  b::Float64 = 2 * tray.x * tray.dx + 2 * tray.y * tray.dy
  c::Float64 = tray.x^2 + tray.y^2 - cylinder.rsq[i]
  d::Float64 = b^2 - (4*a*c)

  t1::Float64 = Inf
  t2::Float64 = Inf
  t::Float64 = Inf

  if d > 0
    t1 = (-1.0*b - sqrt(d))/(2*a)
    t2 = (-1.0*b + sqrt(d))/(2*a)
  end

  # choose the smallest, non-negative solution for t
  if t1 <= 0
    if t2 > 0
      t = t2
    end
  else
    if t2 <= 0
      t = t1
    elseif t1 < t2
      t = t1
    else
      t = t2
    end
  end
  # t1 and t2 can both be negative, hence t can be negative!

  hit::Bool = false
  tI = Point{T}(tray.x + tray.dx * t, tray.y + tray.dy * t, tray.z + tray.dz * t)
  if (t > 0) && (tI.z > 0.0) && (tI.z < cylinder.length[i])
    hit = true
  end

  I = Point{T}(ray.x + ray.dx * t, ray.y + ray.dy * t, ray.z + ray.dz * t)
  nx::Float64, ny::Float64, nz::Float64 = normalize(tI.x, tI.y, 0.0)
  nx, ny, nz = transform(nx, ny, nz, xoff, yoff, zoff, beta, gamma)
  normal::Point = Point{T}(nx, ny, nz)

  return hit, I, normal, t
end


function intersect(mesh::Mesh{T,S,A,E}, i::S, ray::Ray{T}) where {T<:AbstractFloat, S<:Integer, A<:AbstractArray{T,1}, E<:AbstractArray{S,1}}
    w0x = ray.x - mesh.v1x[i]
    w0y = ray.y - mesh.v1y[i]
    w0z = ray.z - mesh.v1z[i]
    a = -((mesh.nx[i] * w0x) + (mesh.ny[i] * w0y) + (mesh.nz[i] * w0z))
    b = (ray.dx * mesh.nx[i]) + (ray.dy * mesh.ny[i]) + (ray.dz * mesh.nz[i])

    if abs(b) < SMALL_NUM
        nullPoint = Point{T}(Inf, Inf, Inf)
        return false, nullPoint, nullPoint, Inf # parallel
    end

    r = a / b
    if r < 0.0
        nullPoint = Point{T}(Inf, Inf, Inf)
        return false, nullPoint, nullPoint, r # away
    end

    # point of plane intersection
    x = ray.x + r * ray.dx
    y = ray.y + r * ray.dy
    z = ray.z + r * ray.dz
    I = Point{T}(x, y, z)
    N = Point{T}(mesh.nx[i], mesh.ny[i], mesh.nz[i])

    mx = x - mesh.v1x[i]
    my = y - mesh.v1y[i]
    mz = z - mesh.v1z[i]
    mu = (mx * mesh.ux[i]) + (my * mesh.uy[i]) + (mz * mesh.uz[i])
    mw = (mx * mesh.wx[i]) + (my * mesh.wy[i]) + (mz * mesh.wz[i])

    s = ((mesh.uw[i] * mw) - (mesh.ww[i] * mu)) / mesh.D[i]
    if (s < 0.0) | (s > 1.0)
        return false, I, N, r # outside
    end

    t = ((mesh.uw[i] * mu) - (mesh.uu[i] * mw)) / mesh.D[i]
    if (t < 0.0) | ((s + t) > 1.0)
        return false, I, N, r # outside
    end

    return true, I, N, r # ray intersects triangle
end


function intersect(geometry::Geometry{T,S,A,E}, ray::Ray{T}) where {T<:AbstractFloat, S<:Integer, A<:AbstractArray{T,1}, E<:AbstractArray{S,1}}
    hit::Bool = false
    In::Point{T} = Point{T}(Inf, Inf, Inf)
    eidx::S = 0
    normal::Point{T} = Point{T}(Inf, Inf, Inf)
    rmin::T = Inf

    # intersect the individual elements of this geometry (e.g. facets)...
    eidxs::Array{S,1} = collect(1:size(geometry.mtl, 1))
    for i::S in eidxs
        hit_, In_, N_, r_ = intersect(geometry, i, ray)
        if (hit_ === true) && (r_ < rmin)
            eidx = i
            rmin = r_
            In = In_
            normal = N_
            hit = true
        end
    end
    return hit, In, eidx, normal, rmin
end

function intersect(geometry::Geometry{T,S,A,E}, bvh::Bvh{T,S,A,E}, ray::Ray{T}) where {T<:AbstractFloat, S<:Integer, A<:AbstractArray{T,1}, E<:AbstractArray{S,1}}
    hit::Bool = false
    In::Point{T} = Point{T}(Inf, Inf, Inf)
    eidx::S = 0
    normal::Point{T} = Point{T}(Inf, Inf, Inf)
    rmin::T = Inf

    # intersect mesh facets
    #eidxs = collect(1:size(geometry.mtl, 1))
    eidxs::Array{S,1} = intersect_bvh(bvh, ray)
    for i::S in eidxs
        # and see if it is nearest
        hit_, I_, N_, r_ = intersect(geometry, i, ray)
        if (hit_ === true) && (r_ < rmin)
            eidx = i
            rmin = r_
            In = I_
            normal = N_
            hit = true
        end
    end
    return hit, In, eidx, normal, rmin
end

function intersect_scene(assets::Array{Asset{T,S},1}, geometries::AbstractArray{Geometry{T,S,A,E},1}, geometry_bvhs::AbstractArray{Bvh{T,S,A,E},1}, scene_bvh::Bvh{T,S,A,E}, ray::Ray{T}) where {T<:AbstractFloat, S<:Integer, A<:AbstractArray{T,1}, E<:AbstractArray{S,1}}
    hit::Bool = false
    In::Point{T} = Point{T}(Inf, Inf, Inf)
    nrml::Point{T} = Point{T}(Inf, Inf, Inf)
    tray::Ray{T} = Ray{T}(Inf, Inf, Inf, Inf, Inf, Inf)
    rmin::T = Inf
    instidx::S = 0
    aidx::S = 0 # asset index
    oidx::S = 0 # object index / geometry index
    eidx::S = 0 # element index (of the geometry indexed with oidx)

    #aidxs::Array{S} = collect(1:length(assets))
    aidxs::Array{S} = intersect_bvh(scene_bvh, ray)
    for i::S in aidxs
        oidx_ = assets[i].oidx

        geometry = geometries[oidx_]
        geometry_bvh = geometry_bvhs[oidx_]

        tray_ = transform_inverse(ray, assets[i]) # from world space to object space
        hit_, In_, eidx_, nrml_, rmin_ = intersect(geometry, geometry_bvh, tray_)
        if (hit_ === true) && (rmin_ < rmin)
            hit = true
            x = ray.x + ray.dx * rmin_
            y = ray.y + ray.dy * rmin_
            z = ray.z + ray.dz * rmin_
            In = Point{T}(x, y, z)
            tray = tray_
            nrml = nrml_
            rmin = rmin_
            aidx = i
            oidx = oidx_
            eidx = eidx_
        end
    end
    return hit, In, tray, nrml, rmin, aidx, oidx, eidx
end
