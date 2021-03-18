function make_sky(light_pos::A, power::Integer, ambient::T, diffusePortion::T) where {T<:AbstractFloat, A<:AbstractArray{T,1}}
  skymap::AbstractArray{T,2} = zeros(3600, 900)
  azi::T = 0
  zen::T = 0
  for i=1:3600
      for j=1:900
          azi = float(i)/1800.0*pi
          zen = float(j)/1800.0*pi
          xi, yi, zi = spherical2carthesian(azi, zen)
          skymap[i,j] = background_radiance([xi, yi, zi], light_pos, power, ambient, diffusePortion)
      end
  end
  return skymap
end

function background_radiance(hemi_sample::Vector{T}, light_pos::Vector{T}, power::Int, ambient::T, diffusePortion) where {T<:AbstractFloat}
  # expand this function with sun position, atmospheric model, etc.
  direct::T = max(0.0, dot(light_pos, hemi_sample))^power
  diffuse::T = max(0.0, dot(light_pos, hemi_sample))^1
  directPortion::T = 1. - diffusePortion
  dirrad::T =  directPortion * direct
  diffrad::T = diffusePortion * diffuse
  radiance::T = (1.0-ambient) * (dirrad + diffrad) + ambient
  return radiance
end

function sample_skymap(skymap::AbstractArray{T, 2}, dx::T, dy::T, dz::T) where {T<:AbstractFloat}
    phi::T, theta::T = carthesian2spherical(dx, dy, dz)
    v::T = phi*1800.0/pi # convert to an index...
    i::Int = round(Int, v)
    w::T = theta*1800.0/pi
    j::Int = round(Int, w)
    if i === 0
        i = 1
    elseif i > 3600
        i = 3600
    end
    if j === 0
        j = 1
    elseif j > 900
        j = 900
    end
    radiance::T = skymap[i, j]
    return radiance
end
