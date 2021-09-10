

function trace_back(geometry::Geometry{T,S,A,E}, shaders::Vector{Function}, ray::Ray{T}, nBands::S, maxNumberOfCycles::S) where {T<:AbstractFloat, S<:Integer, A<:AbstractArray{T,1}, E<:AbstractArray{S,1}}
    # traces ray/photon for rendering
    count_while_loop = 0
    continue_iterating::Bool = true
    radiance::Vector{T} = ones(nBands) .* GREAT_NUM
    mtl::Int = 0
    hit::Bool = false
    while continue_iterating
        count_while_loop += 1
        if count_while_loop >= maxNumberOfCycles
            continue_iterating = false
            continue
        end
        hit_, In, eidx, nrml, rmin = intersect(geometry, ray)
        if hit_ === true
	    hit = true
	    mtl = geometry.mtl[eidx]
            scatter, ray, r = shaders[mtl](In.x, In.y, In.z, -ray.dx, -ray.dy, -ray.dz, nrml.x, nrml.y, nrml.z)
            radiance .*= r
	    if scatter === false
		continue_iterating = false
	    end
        else
	    radiance .*= 0.0
            continue_iterating = false
        end
    end
    return hit, radiance
end

function trace_back(geometry::Geometry{T,S,A,E}, geometry_bvh::Bvh{T,S,A,E}, shaders::Vector{Function}, ray::Ray{T}, nBands::S, maxNumberOfCycles::S) where {T<:AbstractFloat, S<:Integer, A<:AbstractArray{T,1}, E<:AbstractArray{S,1}}
    # traces ray/photon for rendering
    count_while_loop = 0
    continue_iterating::Bool = true
    radiance::Vector{T} = ones(nBands) .* GREAT_NUM
    mtl::Int = 0
    hit::Bool = false
    while continue_iterating
        count_while_loop += 1
        if count_while_loop >= maxNumberOfCycles
            continue_iterating = false
            continue
        end
        hit_, In, eidx, nrml, rmin = intersect(geometry, geometry_bvh, ray)
        if hit_ === true
	    hit = true
	    mtl = geometry.mtl[eidx]
            scatter, ray, r = shaders[mtl](In.x, In.y, In.z, -ray.dx, -ray.dy, -ray.dz, nrml.x, nrml.y, nrml.z)
            radiance .*= r
	    if scatter === false
		continue_iterating = false
	    end
        else
	    radiance .*= 0.0
            continue_iterating = false
        end
    end
    return hit, radiance
end

function trace_back(geometry::Geometry{T,S,A,E}, shaders::Vector{Function}, skymap::B, ray::Ray{T}, nBands::S, maxNumberOfCycles::S) where {T<:AbstractFloat, S<:Integer, A<:AbstractArray{T,1}, B<:AbstractArray{T,2}, E<:AbstractArray{S,1}}
    # traces ray/photon for rendering
    count_while_loop = 0
    continue_iterating::Bool = true
    radiance::Vector{T} = ones(nBands) .* GREAT_NUM
    mtl::Int = 0
    hit::Bool = false
    while continue_iterating
        count_while_loop += 1
        if count_while_loop >= maxNumberOfCycles
            continue_iterating = false
            continue
        end
        hit_, In, eidx, nrml, rmin = intersect(mesh, ray)
        if hit_mesh === true
	    hit = true
	    mtl = geometry.mtl[eidx]
            scatter, ray, r = shaders[mtl](In.x, In.y, In.z, -ray.dx, -ray.dy, -ray.dz, nrml.x, nrml.y, nrml.z)
            radiance .*= r
	    if scatter === false
		continue_iterating = false
	    end
        else
	    radiance .*= sample_skymap(skymap, ray.dx, ray.dy, ray.dz)
            continue_iterating = false
        end
    end
    return hit, radiance
end

function trace_back(geometry::Geometry{T,S,A,E}, geometry_bvh::Bvh{T,S,A,E}, shaders::Vector{Function}, skymap::B, ray::Ray{T}, nBands::S, maxNumberOfCycles::S) where {T<:AbstractFloat, S<:Integer, A<:AbstractArray{T,1}, B<:AbstractArray{T,2}, E<:AbstractArray{S,1}}
    # traces ray/photon for rendering
    count_while_loop = 0
    continue_iterating::Bool = true
    radiance::Vector{T} = ones(nBands) .* GREAT_NUM
    mtl::Int = 0
    hit::Bool = false
    while continue_iterating
        count_while_loop += 1
        if count_while_loop >= maxNumberOfCycles
            continue_iterating = false
            continue
        end
        hit_, In, eidx, nrml, rmin = intersect(geometry, geometry_bvh, ray)
        if hit_ === true
	    hit = true
	    mtl = geometry.mtl[eidx]
            scatter, ray, r = shaders[mtl](In.x, In.y, In.z, -ray.dx, -ray.dy, -ray.dz, nrml.x, nrml.y, nrml.z)
            radiance .*= r
	    if scatter === false
		continue_iterating = false
	    end
        else
	    radiance .*= sample_skymap(skymap, ray.dx, ray.dy, ray.dz)
            continue_iterating = false
        end
    end
    return hit, radiance
end

function trace_back(assets::AbstractArray{Asset{T,S},1}, geometries::AbstractArray{Geometry{T,S,A,E},1}, palettes::AbstractArray{Vector{Function},1}, geometry_bvhs::AbstractArray{Bvh{T,S,A,E},1}, scene_bvh::Bvh{T,S,A,E}, ray::Ray{T}, skymap::B, nBands::S, maxNumberOfCycles::S) where {T<:AbstractFloat, S<:Integer, A<:AbstractArray{T,1}, B<:AbstractArray{T,2}, E<:AbstractArray{S,1}}
    # traces ray/photon for rendering
    count_while_loop::S = 0
    continue_iterating::Bool = true
    radiance::Vector{T} = ones(nBands)
    hit::Bool = false
    while continue_iterating
        count_while_loop += 1
        if count_while_loop >= maxNumberOfCycles
            continue_iterating = false
            continue
        end
        hit_scene, I_scene, tray, nrml, range_scene, aidx, oidx, eidx = intersect_scene(assets, geometries, geometry_bvhs, scene_bvh, ray)
        if hit_scene === true
	    hit = true
            mtl = geometries[oidx].mtl[eidx]
            pidx = assets[aidx].pidx
            shader = palettes[pidx][mtl]
            scatter, ray, r = shader(I_scene.x, I_scene.y, I_scene.z, -tray.dx, -tray.dy, -tray.dz, nrml.x, nrml.y, nrml.z)
            radiance .*= r
	    if scatter === false
		continue_iterating = false
	    end
        else
	    radiance = radiance .* sample_skymap(skymap, ray.dx, ray.dy, ray.dz)
            continue_iterating = false
        end
    end
    return radiance
end



# following traces in forward direction with cyclic Walls. 
# These functions record the "faterecs" and don't use skymaps.
# These function can be used to compute e.g. fAPAR or surface BDRF/BRF/BHR etc.
# but some are relatively slow because of the instantiation of small arrays...

function record(geometry::Geometry{T,S,A,E}, shaders::Vector{Function}, walls::Walls{T}, ray::Ray{T}, maxNumberOfBounces::Int) where {T<:AbstractFloat, S<:Integer, A<:AbstractArray{T,1}, E<:AbstractArray{S,1}}
    # traces ray/photon until it is absorbed or bounces through the reference plane
    #
    # TERMINAL CODES:
    #
    # -1    the maximum number of cycles was reached
    #  1    photon escapes through the reference plane
    #  2    photon escapes through the bottom plane
    #  3    photon cyclicly re-enters the scene (through walls)'
    #  4    photon is absorbed
    #  5    photon is scattered
    #
    path::AbstractArray{T, 2} = zeros(Float64, maxNumberOfBounces+1, 3)
    path[1,:] = [ray.x, ray.y, ray.z]
    codes::AbstractArray{Int, 1} = zeros(Int, maxNumberOfBounces+1) # code[1] === 0
    index::Int = 1
    count_while_loop::Int = 0
    continue_iterating::Bool = true
    mtl::Int = -1
    while continue_iterating
        count_while_loop += 1
        index += 1
        if count_while_loop >= maxNumberOfBounces
            path[index, :] = [ray.x, ray.y, ray.z]
            codes[index] = -1
            continue_iterating = false
            continue
        end
        I_wall, side, r_wall = intersect_walls(walls, ray)
        hit_, In, eidx, nrml, rmin = intersect(geometry, ray)
        wall = r_wall <= rmin
        if wall === true
            if side === 't'
                # photon escapes through reference plane
                continue_iterating = false
                codes[index] = 1
                path[index, :] = [I_wall.x, I_wall.y, I_wall.z]
            elseif (side === 'b')
                # photon is absorbed by a black surface underneath the scene
                continue_iterating = false
                codes[index] = 2
                path[index, :] = [I_wall.x, I_wall.y, I_wall.z]
            else # thus the condition: ' !(side ==='t') && !(side === 'b') ' is true
                # photon cyclicly re-enters from the opposite wall
                ray = cycle_ray(walls, ray)
                count_while_loop -= 1 # don't count these events toward the number of cycles
                index -= 1
            end
        else # implies hit_ is true, because the ray either hits the scene or a wall
            hit = true
            mtl = geometry.mtl[eidx]
            scatter, ray, r = shaders[mtl](In.x, In.y, In.z, -ray.dx, -ray.dy, -ray.dz, nrml.x, nrml.y, nrml.z)
            radiance .*= r
            if scatter === false
                continue_iterating = false
                codes[index] = 4
                path[index, :] = [I_scene.x, I_scene.y, I_scene.z]
            else
                codes[index] = 5
                path[index, :] = [I_scene.x, I_scene.y, I_scene.z]
            end
	end
    end
    return path, codes
end

function record(geometry::Geometry{T,S,A,E}, geometry_bvh::Bvh{T,S,A,E}, shaders::Vector{Function}, walls::Walls{T}, ray::Ray{T}, maxNumberOfBounces::S) where {T<:AbstractFloat, S<:Integer, A<:AbstractArray{T,1}, E<:AbstractArray{S,1}}
    # traces ray/photon until it is absorbed or bounces through the reference plane
    #
    # TERMINAL CODES:
    #
    # -1    the maximum number of cycles was reached
    #  1    photon escapes through the reference plane
    #  2    photon escapes through the bottom plane
    #  3    photon cyclicly re-enters the scene (through walls)'
    #  4    photon is absorbed
    #  5    photon is scattered
    #
    path::AbstractArray{T, 2} = zeros(Float64, maxNumberOfBounces+1, 3)
    path[1,:] = [ray.x, ray.y, ray.z]
    codes::AbstractArray{Int, 1} = zeros(Int, maxNumberOfBounces+1) # code[1] === 0
    index::Int = 1
    count_while_loop::Int = 0
    continue_iterating::Bool = true
    mtl::Int = -1
    while continue_iterating
        count_while_loop += 1
        index += 1
        if count_while_loop >= maxNumberOfBounces
            path[index, :] = [ray.x, ray.y, ray.z]
            codes[index] = -1
            continue_iterating = false
            continue
        end
        I_wall, side, r_wall = intersect_walls(walls, ray)
        hit_, In, eidx, nrml, rmin = intersect(geometry, geometry_bvh, ray)
        wall = r_wall <= rmin
        if wall === true
            if side === 't'
                # photon escapes through reference plane
                continue_iterating = false
                codes[index] = 1
                path[index, :] = [I_wall.x, I_wall.y, I_wall.z]
            elseif (side === 'b')
                # photon is absorbed by a black surface underneath the scene
                continue_iterating = false
                codes[index] = 2
                path[index, :] = [I_wall.x, I_wall.y, I_wall.z]
            else # thus the condition: ' !(side ==='t') && !(side === 'b') ' is true
                # photon cyclicly re-enters from the opposite wall
                ray = cycle_ray(walls, ray)
                count_while_loop -= 1 # don't count these events toward the number of cycles
                index -= 1
            end
        else # implies hit_ is true, because the ray either hits the scene or a wall
            hit = true
            mtl = geometry.mtl[eidx]
            scatter, ray, r = shaders[mtl](In.x, In.y, In.z, -ray.dx, -ray.dy, -ray.dz, nrml.x, nrml.y, nrml.z)
            radiance .*= r
            if scatter === false
                continue_iterating = false
                codes[index] = 4
                path[index, :] = [I_scene.x, I_scene.y, I_scene.z]
            else
                codes[index] = 5
                path[index, :] = [I_scene.x, I_scene.y, I_scene.z]
            end
	end
    end
    return path, codes
end

function record(assets::AbstractArray{Asset{T,S},1}, geometries::AbstractArray{Geometry{T,S,A,E},1}, palettes::AbstractArray{Vector{Function},1}, geometry_bvhs::AbstractArray{Bvh{T,S,A,E},1}, scene_bvh::Bvh{T,S,A,E}, walls::Walls{T}, ray::Ray{T}, maxNumberOfBounces::S) where {T<:AbstractFloat, S<:Integer, A<:AbstractArray{T,1}, E<:AbstractArray{S,1}}
    # traces ray/photon until it is absorbed or bounces through the reference plane
    #
    # TERMINAL CODES:
    #
    # -1    the maximum number of cycles was reached
    #  1    photon escapes through the reference plane
    #  2    photon escapes through the bottom plane
    #  3    photon cyclicly re-enters the scene (through walls)'
    #  4    photon is absorbed
    #  5    photon is scattered
    #
    path::AbstractArray{T, 2} = zeros(Float64, maxNumberOfBounces+1, 3)
    path[1,:] = [ray.x, ray.y, ray.z]
    codes::AbstractArray{Int, 1} = zeros(Int, maxNumberOfBounces+1) # code[1] === 0
    index::Int = 1
    count_while_loop::Int = 0
    continue_iterating::Bool = true
    mtl::Int = -1
    while continue_iterating
        count_while_loop += 1
	index += 1
        if count_while_loop >= maxNumberOfBounces
            path[index, :] = [ray.x, ray.y, ray.z]
            codes[index] = -1
            continue_iterating = false
            continue
        end
        I_wall, side, r_wall = intersect_walls(walls, ray)
        hit_scene, I_scene, tray, nrml, range_scene, aidx, oidx, eidx = intersect_scene(assets, geometries, geometry_bvhs, scene_bvh, ray)
        wall = r_wall <= range_scene
        if wall === true
            if side === 't'
                # photon escapes through reference plane
                continue_iterating = false
                codes[index] = 1
                path[index, :] = [I_wall.x, I_wall.y, I_wall.z]
            elseif (side === 'b')
                # photon is absorbed by a black surface underneath the scene
                continue_iterating = false
                codes[index] = 2
                path[index, :] = [I_wall.x, I_wall.y, I_wall.z]
            else # thus the condition: ' !(side ==='t') && !(side === 'b') ' is true
                # photon cyclicly re-enters from the opposite wall
                ray = cycle_ray(walls, ray)
                count_while_loop -= 1 # don't count these events toward the number of cycles
                index -= 1
            end
        else # wall is false, hence ray must have hit the (objects in the) scene
            hit = true
            @assert(hit_scene)
	    geometry = geometries[oidx]
            mtl = geometry.mtl[eidx]
            pidx = assets[aidx].pidx
            shader = palettes[pidx][mtl]
            scatter, ray, r = shader(I_scene.x, I_scene.y, I_scene.z, -tray.dx, -tray.dy, -tray.dz, nrml.x, nrml.y, nrml.z)
            if scatter === false
                continue_iterating = false
                codes[index] = 4
                path[index, :] = [I_scene.x, I_scene.y, I_scene.z]
            else
                codes[index] = 5
                path[index, :] = [I_scene.x, I_scene.y, I_scene.z]
            end
	end
    end
    return path, codes
end

function trace_forward(geometry::Geometry{T,S,A,E}, shaders::Vector{Function}, walls::Walls{T}, ray::Ray{T}, maxNumberOfBounces::S) where {T<:AbstractFloat, S<:Integer, A<:AbstractArray{T,1}, E<:AbstractArray{S,1}}
    # traces ray/photon until it is absorbed or bounces through the reference plane
    count_while_loop::Int = 0
    continue_iterating::Bool = true
    mtl::Int = -1
    while continue_iterating
        count_while_loop += 1
        if count_while_loop >= maxNumberOfBounces
            continue_iterating = false
            continue
        end
        I_wall, side, r_wall = intersect_walls(walls, ray)
        hit_, In, eidx, nrml, rmin = intersect(geometry, ray)
        wall = r_wall <= rmin
        if wall === true
            if side === 't'
                # photon escapes through reference plane
                continue_iterating = false
            elseif (side === 'b')
                # photon is absorbed by a black surface underneath the scene
                continue_iterating = false
            else # thus the condition: ' !(side ==='t') && !(side === 'b') ' is true
                # photon cyclicly re-enters from the opposite wall
                ray = cycle_ray(walls, ray)
                count_while_loop -= 1 # don't count these events toward the number of cycles
            end
        else # implies hit_ is true, because the ray either hits the scene or a wall
            hit = true
            mtl = geometry.mtl[eidx]
            scatter, ray, r = shaders[mtl](In.x, In.y, In.z, -ray.dx, -ray.dy, -ray.dz, nrml.x, nrml.y, nrml.z)
            radiance .*= r
            if scatter === false
                continue_iterating = false
            end
	end
    end
    return ray
end

function trace_forward(geometry::Geometry{T,S,A,E}, geometry_bvh::Bvh{T,S,A,E}, shaders::Vector{Function}, walls::Walls{T}, ray::Ray{T}, maxNumberOfBounces::S) where {T<:AbstractFloat, S<:Integer, A<:AbstractArray{T,1}, E<:AbstractArray{S,1}}
    # traces ray/photon until it is absorbed or bounces through the reference plane
    count_while_loop::Int = 0
    continue_iterating::Bool = true
    mtl::Int = -1
    while continue_iterating
        count_while_loop += 1
        if count_while_loop >= maxNumberOfBounces
            continue_iterating = false
            continue
        end
        I_wall, side, r_wall = intersect_walls(walls, ray)
        hit_, In, eidx, nrml, rmin = intersect(geometry, geometry_bvh, ray)
        wall = r_wall <= rmin
        if wall === true
            if side === 't'
                # photon escapes through reference plane
                continue_iterating = false
            elseif (side === 'b')
                # photon is absorbed by a black surface underneath the scene
                continue_iterating = false
            else # thus the condition: ' !(side ==='t') && !(side === 'b') ' is true
                # photon cyclicly re-enters from the opposite wall
                ray = cycle_ray(walls, ray)
                count_while_loop -= 1 # don't count these events toward the number of cycles
            end
        else # implies hit_ is true, because the ray either hits the scene or a wall
            hit = true
            mtl = geometry.mtl[eidx]
            scatter, ray, r = shaders[mtl](In.x, In.y, In.z, -ray.dx, -ray.dy, -ray.dz, nrml.x, nrml.y, nrml.z)
            radiance .*= r
            if scatter === false
                continue_iterating = false
            end
	end
    end
    return ray
end

function trace_forward(assets::AbstractArray{Asset{T,S},1}, geometries::AbstractArray{Geometry{T,S,A,E},1}, palettes::AbstractArray{Vector{Function},1}, geometry_bvhs::AbstractArray{Bvh{T,S,A,E},1}, scene_bvh::Bvh{T,S,A,E}, walls::Walls{T}, ray::Ray{T}, maxNumberOfBounces::S) where {T<:AbstractFloat, S<:Integer, A<:AbstractArray{T,1}, E<:AbstractArray{S,1}}
    # traces ray/photon until it is absorbed or bounces through the reference plane
    count_while_loop::Int = 0
    continue_iterating::Bool = true
    mtl::Int = -1
    while continue_iterating
        count_while_loop += 1
        if count_while_loop >= maxNumberOfBounces
            push!(codes, -1)
            push!(path, Line(ray.x, ray.y, ray.z, ray.x, ray.y, ray.z)) # add a zero-length line, basically indicating that we stop tracing here...
            continue_iterating = false
            continue
        end
        I_wall, side, r_wall = intersect_walls(walls, ray)
        hit_scene, I_scene, tray, nrml, range_scene, aidx, oidx, eidx = intersect_scene(assets, geometries, geometry_bvhs, scene_bvh, ray)
        wall = r_wall <= range_scene
        if wall === true
            if side === 't'
                # photon escapes through reference plane
                continue_iterating = false
            elseif (side === 'b')
                # photon is absorbed by a black surface underneath the scene
                continue_iterating = false
            else # thus the condition: ' !(side ==='t') && !(side === 'b') ' is true
                # photon cyclicly re-enters from the opposite wall
                ray = cycle_ray(walls, ray)
                count_while_loop -= 1 # don't count these events toward the number of cycles
            end
        else # wall is false, hence ray must have hit the (objects in the) scene
            hit = true
	    geometry = geometries[oidx]
            mtl = geometry.mtl[eidx]
            pidx = assets[aidx].pidx
            shader = palettes[pidx][mtl]
            scatter, ray, r = shader(I_scene.x, I_scene.y, I_scene.z, -tray.dx, -tray.dy, -tray.dz, nrml.x, nrml.y, nrml.z)
            if scatter === false
                continue_iterating = false
            end
	end
    end
    return ray
end
