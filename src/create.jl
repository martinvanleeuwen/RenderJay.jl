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


