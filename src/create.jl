function create_ray(x1::T, y1::T, z1::T, x2::T, y2::T, z2::T) where {T<:AbstractFloat}
    dx = x2 - x1
    dy = y2 - y1
    dz = z2 - z1
    l = sqrt((dx^dx) + (dy*dy) + (dz*dz))
    dx = dx / l
    dy = dy / l
    dz = dz / l
    ray = Ray(x1, y1, z1, dx, dy, dz)
    return ray
end

function create_rays(minx::AbstractFloat, miny::AbstractFloat, maxx::AbstractFloat, maxy::AbstractFloat, height::AbstractFloat, dx::AbstractFloat, dy::AbstractFloat, dz::AbstractFloat, n::Int)
    x = rand(Uniform(minx, maxx), n)
    y = rand(Uniform(miny, maxy), n)
    rays = Ray[]
    for i=1:n
        push!(rays, Ray(x[i], y[i], height, dx, dy, dz))
    end
    return rays
end

function create_ray(walls::Walls, dx::AbstractFloat, dy::AbstractFloat, dz::AbstractFloat)
    minx, maxx = walls.west + SMALL_NUM, walls.east - SMALL_NUM
    miny, maxy = walls.south + SMALL_NUM, walls.north - SMALL_NUM # check this though!
    height = walls.top - SMALL_NUM
    x = rand(Uniform(minx, maxx))
    y = rand(Uniform(miny, maxy))
    ray = Ray(x, y, height, dx, dy, dz)
    return ray
end

function create_rays(walls::Walls, dx::AbstractFloat, dy::AbstractFloat, dz::AbstractFloat, n::Int)
    minx, maxx = walls.west + SMALL_NUM, walls.east - SMALL_NUM
    miny, maxy = walls.south + SMALL_NUM, walls.north - SMALL_NUM # check this though!
    height = walls.top - SMALL_NUM
    rays = create_rays(minx, miny, maxx, maxy, height, dx, dy, dz, n)
    return rays
end

function create_coords(camera::Camera{T,S,A}) where {T<:AbstractFloat, S<:Integer, A<:AbstractArray{T,1}}
    # construct an array of coordinate, just for the shape...
    coords = Array{Coord{S}}(undef, camera.xResolution, camera.yResolution)
    for x::S=1:camera.xResolution
        for y::S=1:camera.yResolution
            coords[x,y] = Coord(x, y)
        end
    end
    return coords
end

