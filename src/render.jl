

function render_pixel(coord::Coord{S}, geometry::Geometry{T,S,A,E}, shaders::Vector{Function}, camera::Camera{T,S,A}) where {T<:AbstractFloat, S<:Integer, A<:AbstractArray{T,1}, E<:AbstractArray{S,1}}
    # check possible speed-ups by preallocating temporary arrays...
    xr::T = 0.0
    xrxIV::Vector{T} = [0.0, 0.0, 0.0]
    yr::T = 0.0
    yryIV::Vector{T} = [0.0, 0.0, 0.0]
    radiance::Vector{T} = zeros(camera.nBands)
    for i = 1:camera.nRaysPerPixel
        xr = coord.x + rand()
        xrxIV = xr .* camera.xIncVector
        yr = coord.y + rand()
        yryIV = yr .* camera.yIncVector
        viewPlanePoint = camera.viewPlaneBottomLeftPoint .+ xrxIV .+ yryIV
        dx, dy, dz = viewPlanePoint .- camera.eyePoint
        ray = Ray(camera.eyePoint[1], camera.eyePoint[2], camera.eyePoint[3], dx, dy, dz)
	hit, r = trace_back(geometry, shaders, ray, camera.nBands, camera.maxNumberOfCycles)
        radiance = radiance .+ r
    end
    n::T = camera.nRaysPerPixel
    radiance = radiance ./ n
    return radiance
end


function render_pixel(coord::Coord{S}, geometry::Geometry{T,S,A,E}, shaders::Vector{Function}, skymap::B,  camera::Camera{T,S,A}) where {T<:AbstractFloat, S<:Integer, A<:AbstractArray{T,1}, B<:AbstractArray{T,2}, E<:AbstractArray{S,1}}
    # check possible speed-ups by preallocating temporary arrays...
    xr::T = 0.0
    xrxIV::Vector{T} = [0.0, 0.0, 0.0]
    yr::T = 0.0
    yryIV::Vector{T} = [0.0, 0.0, 0.0]
    radiance::Vector{T} = zeros(camera.nBands)
    for i = 1:camera.nRaysPerPixel
        xr = coord.x + rand()
        xrxIV = xr .* camera.xIncVector
        yr = coord.y + rand()
        yryIV = yr .* camera.yIncVector
        viewPlanePoint = camera.viewPlaneBottomLeftPoint .+ xrxIV .+ yryIV
        dx, dy, dz = viewPlanePoint .- camera.eyePoint
        ray = Ray(camera.eyePoint[1], camera.eyePoint[2], camera.eyePoint[3], dx, dy, dz)
	hit, r = trace_back(geometry, shaders, skymap, ray, camera.nBands, camera.maxNumberOfCycles)
        radiance = radiance .+ r
    end
    n::T = camera.nRaysPerPixel
    radiance = radiance ./ n
    return radiance
end


function render_pixel(coord::Coord{S}, geometry::Geometry{T,S,A,E}, geometry_bvh::Bvh{T,S,A,E}, shaders::Vector{Function}, camera::Camera{T,S,A}) where {T<:AbstractFloat, S<:Integer, A<:AbstractArray{T,1}, E<:AbstractArray{S,1}}
    # check possible speed-ups by preallocating temporary arrays...
    xr::T = 0.0
    xrxIV::Vector{T} = [0.0, 0.0, 0.0]
    yr::T = 0.0
    yryIV::Vector{T} = [0.0, 0.0, 0.0]
    radiance::Vector{T} = zeros(camera.nBands)
    for i = 1:camera.nRaysPerPixel
        xr = coord.x + rand()
        xrxIV = xr .* camera.xIncVector
        yr = coord.y + rand()
        yryIV = yr .* camera.yIncVector
        viewPlanePoint = camera.viewPlaneBottomLeftPoint .+ xrxIV .+ yryIV
        dx, dy, dz = viewPlanePoint .- camera.eyePoint
        ray = Ray(camera.eyePoint[1], camera.eyePoint[2], camera.eyePoint[3], dx, dy, dz)
	hit, r = trace_back(geometry, geometry_bvh, shaders, ray, camera.nBands, camera.maxNumberOfCycles)
	radiance = radiance .+ r
    end
    n::T = camera.nRaysPerPixel
    radiance = radiance ./ n
    return radiance
end


function render_pixel(coord::Coord{S}, assets::AbstractArray{Asset{T,S},1}, geometries::AbstractArray{Geometry{T,S,A,E},1}, palettes::AbstractArray{Vector{Function},1}, geometry_bvhs::AbstractArray{Bvh{T,S,A,E},1}, scene_bvh::Bvh{T,S,A,E}, skymap::B, camera::Camera{T,S,A}) where {T<:AbstractFloat, S<:Integer, A<:AbstractArray{T,1}, B<:AbstractArray{T,2}, E<:AbstractArray{S,1}}
    # check possible speed-ups by preallocating temporary arrays...
    radiance::Vector{T} = zeros(camera.nBands)
    for i = 1:camera.nRaysPerPixel
        xr = coord.x + rand()
        xrxIVx = xr * camera.xIncVector[1]
        xrxIVy = xr * camera.xIncVector[2]
        xrxIVz = xr * camera.xIncVector[3]
        yr = coord.y + rand()
        yryIVx = yr * camera.yIncVector[1]
        yryIVy = yr * camera.yIncVector[2]
        yryIVz = yr * camera.yIncVector[3]
        viewPlanePointx = camera.viewPlaneBottomLeftPoint[1] + xrxIVx + yryIVx
        viewPlanePointy = camera.viewPlaneBottomLeftPoint[2] + xrxIVy + yryIVy
        viewPlanePointz = camera.viewPlaneBottomLeftPoint[3] + xrxIVz + yryIVz
        dx = viewPlanePointx - camera.eyePoint[1]
	dy = viewPlanePointy - camera.eyePoint[2]
	dz = viewPlanePointz - camera.eyePoint[3]
        ray = Ray(camera.eyePoint[1], camera.eyePoint[2], camera.eyePoint[3], dx, dy, dz)
        radiance .+= trace_back(assets, geometries, palettes, geometry_bvhs, scene_bvh, ray, skymap, camera.nBands, camera.maxNumberOfCycles)
    end
    n::T = camera.nRaysPerPixel
    radiance = radiance ./ n
    return radiance
end


# although these render_image functions can be used fine to produce the right results, pmap
# is not that quick in this case... It is better to use @sync @distributed in combination with
# SharedArrays, like shown in the testRenderingWytham jupyter notebook...

function render_image(geometry::Geometry{T,S,A,E}, shaders::Vector{Function}, camera::Camera{T,S,A}) where {T<:AbstractFloat, S<:Integer, A<:AbstractArray{T,1}, E<:AbstractArray{S,1}}
    coords = create_coords(camera)
    img = @showprogress pmap(coord -> render_pixel(coord, geometry, shaders, camera), coords)
    return img
end

function render_image(geometry::Geometry{T,S,A,E}, shaders::Vector{Function}, skymap::B, camera::Camera{T,S,A}) where {T<:AbstractFloat, S<:Integer, A<:AbstractArray{T,1}, B<:AbstractArray{T,2}, E<:AbstractArray{S,1}}
    coords = create_coords(camera)
    img = @showprogress pmap(coord -> render_pixel(coord, geometry, shaders, skymap, camera), coords)
    return img
end

function render_image(geometry::Geometry{T,S,A,E}, geometry_bvh::Bvh{T,S,A,E}, shaders::Vector{Function}, camera::Camera{T,S,A}) where {T<:AbstractFloat, S<:Integer, A<:AbstractArray{T,1}, E<:AbstractArray{S,1}}
    coords = create_coords(camera)
    img = @showprogress pmap(coord -> render_pixel(coord, geometry, geometry_bvh, shaders, camera), coords)
    return img
end

function render_image(assets::AbstractArray{Asset{T,S},1}, geometries::AbstractArray{Geometry{T,S,A,E},1}, palette::AbstractArray{Vector{Function},1}, geometry_bvhs::AbstractArray{Bvh{T,S,A,E},1}, scene_bvh::Bvh{T,S,A,E}, skymap::B, camera::Camera{T,S,A}) where {T<:AbstractFloat, S<:Integer, A<:AbstractArray{T,1}, B<:AbstractArray{T,2}, E<:AbstractArray{S,1}}
    coords = create_coords(camera)
    img = @showprogress pmap(coord -> render_pixel(coord, assets, geometries, palette, geometry_bvhs, scene_bvh, skymap, camera), coords)
    return img
end

