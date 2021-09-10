
function get_extremes(x::Vector{T}, y::Vector{T}, z::Vector{T}) where {T<:AbstractFloat}
    minx, maxx = minimum(x), maximum(x)
    miny, maxy = minimum(y), maximum(y)
    minz, maxz = minimum(z), maximum(z)
    p::Vector{T} = [minx, minx, minx, minx, maxx, maxx, maxx, maxx]
    q::Vector{T} = [miny, miny, maxy, maxy, miny, miny, maxx, maxx]
    r::Vector{T} = [minz, maxz, minz, maxz, minz, maxz, minz, maxz]
    return p, q, r
end

function get_translation_parameters(hook::Hook{T,S}) where {T<:AbstractFloat, S<:Integer}
    x = hook.xoff
    y = hook.yoff
    z = hook.zoff
    return x, y, z
end

function get_translation_parameters(transform::Transform{T,S}) where {T<:AbstractFloat, S<:Integer}
    x = transform.matrix[1,4]
    y = transform.matrix[2,4]
    z = transform.matrix[3,4]
    return x, y, z
end


function get_extremes(assets::Array{Asset{T,S}, 1}, geometry_bvhs::Array{Bvh{T,S,A,E}, 1}) where {T<:AbstractFloat, S<:Integer, A<:AbstractArray{T,1}, E<:AbstractArray{S,1}}
    n = size(assets, 1)
    boxes::Array{Float64, 2} = zeros(6, n)
    for i=1:n
        p = assets[i].oidx
	
        xmin = geometry_bvhs[p].xmin[1]
        ymin = geometry_bvhs[p].ymin[1]
        zmin = geometry_bvhs[p].zmin[1]
	
        xmax = geometry_bvhs[p].xmax[1]
        ymax = geometry_bvhs[p].ymax[1]
        zmax = geometry_bvhs[p].zmax[1]
	
        x1, y1, z1 = transform_forward(xmin, ymin, zmin, assets[i])
        x2, y2, z2 = transform_forward(xmin, ymin, zmax, assets[i])
        x3, y3, z3 = transform_forward(xmin, ymax, zmin, assets[i])
        x4, y4, z4 = transform_forward(xmin, ymax, zmax, assets[i])
        x5, y5, z5 = transform_forward(xmax, ymin, zmin, assets[i])
        x6, y6, z6 = transform_forward(xmax, ymin, zmax, assets[i])
        x7, y7, z7 = transform_forward(xmax, ymax, zmin, assets[i])
        x8, y8, z8 = transform_forward(xmax, ymax, zmax, assets[i])
	
        txmin = min(x1, x2, x3, x4, x5, x6, x7, x8)
        tymin = min(y1, y2, y3, y4, y5, y6, y7, y8)
        tzmin = min(z1, z2, z3, z4, z5, z6, z7, z8)
	
        txmax = max(x1, x2, x3, x4, x5, x6, x7, x8)
        tymax = max(y1, y2, y3, y4, y5, y6, y7, y8)
        tzmax = max(z1, z2, z3, z4, z5, z6, z7, z8)
	
        boxes[:,i] = [txmin, tymin, tzmin, txmax, tymax, tzmax]
    end
    xmin, xmax = minimum(boxes[1,:]), maximum(boxes[4,:])
    ymin, ymax = minimum(boxes[2,:]), maximum(boxes[5,:])
    zmin, zmax = minimum(boxes[3,:]), maximum(boxes[6,:])
    x = vcat(xmin, xmax)
    y = vcat(ymin, ymax)
    z = vcat(zmin, zmax)
    p, q, r = get_extremes(x, y, z)
    return p, q, r
end



# these functions below are not particularly accurate, so use with caution! They do not consider diameters of disks, cones, cylinders, or assemblies, and no orientations/scaling...

function get_extremes(disk::Disk{T,S,A,E}) where {T<:AbstractFloat, S<:Integer, A<:AbstractArray{T,1}, E<:AbstractArray{S,1}}
    x = disk.p1x
    y = disk.p1y
    z = disk.p1z
    p, q, r = get_extremes(x, y, z)
    return p, q, r
end

function get_extremes(cone::Cone{T,S,A,E}) where {T<:AbstractFloat, S<:Integer, A<:AbstractArray{T,1}, E<:AbstractArray{S,1}}
    x = cone.m1x
    y = cone.m1y
    z = cone.m1z
    p, q, r = get_extremes(x, y, z)
    return p, q, r
end

function get_extremes(cylinder::Cylinder{T,S,A,E}) where {T<:AbstractFloat, S<:Integer, A<:AbstractArray{T,1}, E<:AbstractArray{S,1}}
    x = cylinder.m1x
    y = cylinder.m1y
    z = cylinder.m1z
    p, q, r = get_extremes(x, y, z)
    return p, q, r
end

function get_extremes(mesh::Mesh{T,S,A,E}) where {T<:AbstractFloat, S<:Integer, A<:AbstractArray{T,1}, E<:AbstractArray{S,1}}
    x = vcat(mesh.v1x, mesh.v2x, mesh.v3x)
    y = vcat(mesh.v1y, mesh.v2y, mesh.v3y)
    z = vcat(mesh.v1z, mesh.v2z, mesh.v3z)
    p, q, r = get_extremes(x, y, z)
    return p, q, r
end

function get_extremes(geometry::Geometry{T,S,A,E}, hook::Hook{T,S}) where {T<:AbstractFloat, S<:Integer, A<:AbstractArray{T,1}, E<:AbstractArray{S,1}}
    x, y, z = get_extremes(geometry)
    x, y, z = transform_forward(x, y, z, hook)
    x, y, z = get_extremes(x, y, z)
    return x, y, z
end









#function get_extremes(assets::Array{Asset{T,S},1}, geometry_bvhs::AbstractArray{Bvh{T,S,A,E},1}) where {T<:AbstractFloat, S<:Integer, A<:AbstractArray{T,1}, E<:AbstractArray{S,1}}
#    n = size(assets,1)
#    x = Float64[]
#    y = Float64[]
#    z = Float64[]
#    for i=1:n
#        asset = assets[i]
#        bvh = geometry_bvhs[asset.oidx]
#        push!(x, bvh.xmin...)
#        push!(x, bvh.xmax...)
#        push!(y, bvh.ymin...)
#        push!(y, bvh.ymax...)
#        push!(z, bvh.zmin...)
#        push!(z, bvh.zmax...)
#    end
#    p, q, r = get_extremes(x, y, z)
#    return p, q, r
#end

#function get_extremes(assets::Array{Asset{T,S},1}) where {T<:AbstractFloat, S<:Integer}
#    n = size(assets,1)
#    x = Float64[]
#    y = Float64[]
#    z = Float64[]
#    for i=1:n
#        asset = assets[i]
#        a, b, c = get_translation_parameters(asset)
#        push!(x, a)
#        push!(y, b)
#        push!(z, c)
#    end
#    p, q, r = get_extremes(x, y, z)
#    return p, q, r
#end
