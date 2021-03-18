

function make_boxes(box::Vector{T}, n::S, scale_factor::T) where {T<:AbstractFloat, S<:Integer}
    # Generates data for testing bvh, returns a set of n boxes inside box
    minx, miny, minz, maxx, maxy, maxz = box
    dx::T = abs(maxx - minx) * (1.0 - scale_factor)
    dy::T = abs(maxy - miny) * (1.0 - scale_factor)
    dz::T = abs(maxz - minz) * (1.0 - scale_factor)
    sx::T = abs(maxx - minx) * scale_factor
    sy::T = abs(maxy - miny) * scale_factor
    sz::T = abs(maxz - minz) * scale_factor
    newBoxes::Array{T,2} = zeros(6, n)
    for i=1:n
        a = rand() * dx + minx
        b = rand() * dy + miny
        c = rand() * dz + minz
        newBoxes[:,i] = [a, b, c, a+sx, b+sy, c+sz]
    end
    return newBoxes
end


function produce_bvh_data(bounding_boxes::Array{T,2}) where T<:AbstractFloat
  # BVH with one leaf per end node
  n::Int = size(bounding_boxes, 2)
  m::Int = 2*n-1
  outerbox::Vector{T} = build_containing_box(bounding_boxes)
  bvh_boxes::Array{T,2} = zeros(6, m)
  bvh_boxes[:,1] = outerbox
  bb2bvh::Array{Int,1} = ones(n) # bounding boxes pointing to their node in the bvh
  bvh2bb::Array{Int,1} = zeros(m)
  bb_idx::Array{Int,1} = collect(1:n)
  left_child::Array{Int} = ones(m)
  right_child::Array{Int} = ones(m)

  writeCursor::Int = 1
  #axis::Char = 'x'
  cutDirection::Array{Char} = ['x' for i::Int = 1:m]
  @showprogress for readCursor::Int = 1:m # sequence through the BVH nodes top-down and slice, if multiple leaves are left within
    #if readCursor % 1000 === 0
    #    println(readCursor, " / ", m)
    #end

    subset_idx = bb_idx[bb2bvh .== readCursor]

    if size(subset_idx, 1) == 0
        continue
    end

    if size(subset_idx, 1) == 1
        bvh2bb[readCursor] = subset_idx[1]
    end

    if size(subset_idx, 1) > 1
      axis = cutDirection[readCursor]
      idxs_left, idxs_right, axis = median_cut(bounding_boxes[:,subset_idx]; axis=axis)
      idxs_left, idxs_right = subset_idx[idxs_left], subset_idx[idxs_right]

      writeCursor += 1
      bvh_boxes[:,writeCursor] = build_containing_box(bounding_boxes[:,idxs_left])
      bb2bvh[idxs_left] .= writeCursor
      left_child[readCursor] = writeCursor
      cutDirection[writeCursor] = axis

      writeCursor += 1
      bvh_boxes[:,writeCursor] = build_containing_box(bounding_boxes[:,idxs_right])
      bb2bvh[idxs_right] .= writeCursor
      right_child[readCursor] = writeCursor
      cutDirection[writeCursor] = axis
    end
  end
  # 1D arrays seem to be faster than a 2D array, even if column-major...
  minx::Array{T,1} = bvh_boxes[1,:]
  miny::Array{T,1} = bvh_boxes[2,:]
  minz::Array{T,1} = bvh_boxes[3,:]
  maxx::Array{T,1} = bvh_boxes[4,:]
  maxy::Array{T,1} = bvh_boxes[5,:]
  maxz::Array{T,1} = bvh_boxes[6,:]
  return minx, miny, minz, maxx, maxy, maxz, bvh2bb, left_child, right_child, bb2bvh
end


function median_cut(bounding_boxes::Array{T,2}; axis::Char=x) where {T<:AbstractFloat}
    n::Int = size(bounding_boxes,2)
    a::Vector{T} = zeros(n)
    if axis == 'x'
        for i::Int = 1:n
            a[i] = bounding_boxes[1,i]
        end
        s = sortperm(a)
        m = floor(Int, size(a,1) / 2)
        idx_left = s[1:m]
        idx_right = s[m+1:end]
        axis = 'y'
    elseif axis == 'y'
        for i::Int = 1:n
            a[i] = bounding_boxes[2,i]
        end
        s = sortperm(a)
        m = floor(Int, size(a,1) / 2)
        idx_left = s[1:m]
        idx_right = s[m+1:end]
        axis = 'z'
    elseif axis == 'z'
        for i::Int = 1:n
            a[i] = bounding_boxes[3,i]
        end
        s = sortperm(a)
        m = floor(Int, size(a,1) / 2)
        idx_left = s[1:m]
        idx_right = s[m+1:end]
        axis = 'x'
    end
    return idx_left, idx_right, axis
end


function build_containing_box(boxes::Array{T,2}) where {T<:AbstractFloat}
    minx::T = Inf
    miny::T = Inf
    minz::T = Inf
    maxx::T = -1.0*Inf
    maxy::T = -1.0*Inf
    maxz::T = -1.0*Inf
    n::Int = size(boxes,2)

    for i::Int = 1:n
        if boxes[1,i] < minx # minPoint.x
            minx = boxes[1,i]
        end
        if boxes[4,i] > maxx # maxPoint.x
            maxx = boxes[4,i]
        end
        if boxes[2,i] < miny # minPoint.y
            miny = boxes[2,i]
        end
        if boxes[5,i] > maxy # maxPoint.y
            maxy = boxes[5,i]
        end
        if boxes[3,i] < minz # minPoint.z
            minz = boxes[3,i]
        end
        if boxes[6,i] > maxz #  maxPoint.z
            maxz = boxes[6,i]
        end
    end
    box::Vector{T} = [minx, miny, minz, maxx, maxy, maxz]
    return box
end


function produce_bounding_boxes(mesh::Mesh)
    n = size(mesh.v1x, 1)
    boxes::Array{Float32, 2} = zeros(6, n)
    for i=1:n
        minx = min(mesh.v1x[i], mesh.v2x[i], mesh.v3x[i])
        miny = min(mesh.v1y[i], mesh.v2y[i], mesh.v3y[i])
        minz = min(mesh.v1z[i], mesh.v2z[i], mesh.v3z[i])
        maxx = max(mesh.v1x[i], mesh.v2x[i], mesh.v3x[i])
        maxy = max(mesh.v1y[i], mesh.v2y[i], mesh.v3y[i])
        maxz = max(mesh.v1z[i], mesh.v2z[i], mesh.v3z[i])
        boxes[:,i] = [minx, miny, minz, maxx, maxy, maxz]
    end
    return boxes
end


function produce_bounding_boxes(cylinder::Cylinder)
    n = size(cylinder.m1x, 1)
    boxes::Array{Float32, 2} = zeros(6, n)
    for i=1:n
        minx = min(cylinder.m1x[i], cylinder.m2x[i])
        miny = min(cylinder.m1y[i], cylinder.m2y[i])
        minz = min(cylinder.m1z[i], cylinder.m2z[i])
        maxx = max(cylinder.m1x[i], cylinder.m2x[i])
        maxy = max(cylinder.m1y[i], cylinder.m2y[i])
        maxz = max(cylinder.m1z[i], cylinder.m2z[i])
        r = cylinder.radius[i]
        boxes[:,i] = [minx-r, miny-r, minz-r, maxx+r, maxy+r, maxz+r]
    end
    return boxes
end


function produce_bounding_boxes(cone::Cone)
    n = size(cone.m1x, 1)
    boxes::Array{Float32, 2} = zeros(6, n)
    for i=1:n
        minx = min(cone.m1x[i] - cone.radius_base[i], cone.v2x[i] - cone.radius_tip[i])
        miny = min(cone.m1y[i] - cone.radius_base[i], cone.v2y[i] - cone.radius_tip[i])
        minz = min(cone.m1z[i] - cone.radius_base[i], cone.v2z[i] - cone.radius_tip[i])
        maxx = max(cone.m1x[i] + cone.radius_base[i], cone.v2x[i] + cone.radius_tip[i])
        maxy = max(cone.m1y[i] + cone.radius_base[i], cone.v2y[i] + cone.radius_tip[i])
        maxz = max(cone.m1z[i] + cone.radius_base[i], cone.v2z[i] + cone.radius_tip[i])
        boxes[:,i] = [minx, miny, minz, maxx, maxy, maxz]
    end
    return boxes
end


function make_bvh(bb::AbstractArray{T}) where {T<:AbstractFloat}
    # create bvh from bounding boxes (6xN array)
    minx, miny, minz, maxx, maxy, maxz, bvh2bb, left_child, right_child, bb2bvh = produce_bvh_data(bb)
    bvh = Bvh(minx, miny, minz, maxx, maxy, maxz, bvh2bb, left_child, right_child)
    return bvh
end


function make_bvh(geometry::Geometry{T,S,A,E}) where {T<:AbstractFloat, S<:Integer, A<:AbstractArray, E<:AbstractArray}
    # use this function only for meshes; scenes are handled below...
    boxes = produce_bounding_boxes(geometry)
    minx, miny, minz, maxx, maxy, maxz, bvh2bb, left_child, right_child, bb2bvh = produce_bvh_data(boxes)
    bvh = Bvh(minx, miny, minz, maxx, maxy, maxz, bvh2bb, left_child, right_child)
    return bvh
end


function make_bvh(assets::Array{Asset{Float64, Int64},1}, mesh_bvhs::Array{Bvh{Float64, Int64, Array{Float64,1}, Array{Int64,1}},1})
    n = size(assets, 1)
    boxes::Array{Float64, 2} = zeros(6, n)
    for i=1:n
        p = assets[i].oidx
        xoff = assets[i].xoff
        yoff = assets[i].yoff
        zoff = assets[i].zoff
        beta = assets[i].beta
        gamma = assets[i].gamma

        xmin = mesh_bvhs[p].xmin[1]
        ymin = mesh_bvhs[p].ymin[1]
        zmin = mesh_bvhs[p].zmin[1]

        xmax = mesh_bvhs[p].xmax[1]
        ymax = mesh_bvhs[p].ymax[1]
        zmax = mesh_bvhs[p].zmax[1]

        x1, y1, z1 = transform(xmin, ymin, zmin, xoff, yoff, zoff, beta, gamma)
        x2, y2, z2 = transform(xmin, ymin, zmax, xoff, yoff, zoff, beta, gamma)
        x3, y3, z3 = transform(xmin, ymax, zmin, xoff, yoff, zoff, beta, gamma)
        x4, y4, z4 = transform(xmin, ymax, zmax, xoff, yoff, zoff, beta, gamma)
        x5, y5, z5 = transform(xmax, ymin, zmin, xoff, yoff, zoff, beta, gamma)
        x6, y6, z6 = transform(xmax, ymin, zmax, xoff, yoff, zoff, beta, gamma)
        x7, y7, z7 = transform(xmax, ymax, zmin, xoff, yoff, zoff, beta, gamma)
        x8, y8, z8 = transform(xmax, ymax, zmax, xoff, yoff, zoff, beta, gamma)

        txmin = min(x1, x2, x3, x4, x5, x6, x7, x8)
        tymin = min(y1, y2, y3, y4, y5, y6, y7, y8)
        tzmin = min(z1, z2, z3, z4, z5, z6, z7, z8)

        txmax = max(x1, x2, x3, x4, x5, x6, x7, x8)
        tymax = max(y1, y2, y3, y4, y5, y6, y7, y8)
        tzmax = max(z1, z2, z3, z4, z5, z6, z7, z8)

        boxes[:,i] = [txmin, tymin, tzmin, txmax, tymax, tzmax]
    end
    minx, miny, minz, maxx, maxy, maxz, bvh2bb, left_child, right_child, bb2bvh = produce_bvh_data(boxes)
    scene_bvh::Bvh{Float64, Int64, Array{Float64,1}, Array{Int64,1}} = Bvh(minx, miny, minz, maxx, maxy, maxz, bvh2bb, left_child, right_child)
    return scene_bvh
end


