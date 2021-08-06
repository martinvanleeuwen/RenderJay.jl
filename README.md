# RenderJay.jl
a Julia-based path tracer

RenderJay is a Julia-based path tracer that is intended for operation on large compute clusters, with applications in ecology and remote sensing where scenes are characterized by very high polygon counts, lots of detail, and light spectra may be broken down into hundreds of wavebands (imaging spectroscopy), but not so much for applications in computer graphics. For example, RenderJay uses no bump maps or textures; everything is down to geometric shaped (i.e., spheres, cylinders, cones, disks and triangle meshes) and bidirectional reflectance (transmittance) distribution functions (BRDFs, BTDFs). Being written in Julia it can profit from Julia's distributed computing capabilities, enabling processing of large workloads across multiple servers.

RenderJay can be installed as follows from the Julia REPL:

```
using Pkg
Pkg.add("RenderJay")
```

# Usage

The following is an example code for rendering. It will produce a top-down view of four Cornell boxes in different colours, floating above a gray flat surface. Already produced renderings can be found in the img/ folder. The surface geometry, downwelling irradiance, as well as the scene specification (XML) file with the shaders can all be found in the test/ folder.

The Wytham Woods image is an example rendering that was derived from data that was provided through the RAdiative transfer Model Intercomparison (RAMI) phase-V, that can be found here:

https://rami-benchmark.jrc.ec.europa.eu/_www/phase_descr.php?strPhase=RAMI5
(feel free to contact Martin van Leeuwen for conversion scripts to Jay format)


![image Wytham Woods](https://github.com/martinvanleeuwen/RenderJay.jl/blob/main/img/wytham.png)


# Plant models

A good place to look for plant models is online at various 3D model warehouses. The Wavefront OBJ file format is very close to what RenderJay uses and a simple script in your favorite language can be used to produce the separate vertex and mesh connectivity files. To see how they should look, check out the assets you find under the test/ folder in this repository. Don't worry about the BVH file (it stands for Bounding Volume Hierarchy and it is for spatial indexing). These are produced the first time you use a new asset. If you ever make changes to an existing asset, do remove the BVH file that is associated with that asset. This will automatically create a new BVH based on the changes you applied to that asset. Besides browsing any 3D warehouses, you can also try out e.g. Arbaro - a tree generator for PovRay that was developed by Weber and Penn 1995. It can be found here: http://arbaro.sourceforge.net/


# Performance indication

Rendering a 512x512px image of the below Cornell boxes scene took 00:34:04 (HH:MM:SS) on the Dell T7910 (DUAL E5-2630V3) using 30 workers and it took 3:39:20 on a Lenovo Edge 15 laptop (i7-4510U) with only 3 workers. For some very simple scenes with only a few cylinders on a flat surface (not shown) up to a million rays per second were processed using 30 workers.


# Application

See this paper (https://doi.org/10.1016/j.rse.2021.112405) for an idea of the kind of cool things that you can do with this module.


# RAMI-5

As part of the RAdiative transfer Model Intercomparison excersize phase 5 we have prepared RenderJay scene files for the various test cases involved, including abstract scenes, fully parameterized canopies, and (semi-)empirical forest canopies. Here is a selection of nadir renderings including some of these scenes. No white-balancing has been applied yet, which causes the Wytham woods scene to appear somewhat greenish. The irradiance conditions for these example renderings were set to a 10% diffuse sky with a little bit of a watery sun.

![image RAMI-5](https://github.com/martinvanleeuwen/RenderJay.jl/blob/main/img/rami5_test.png)

Here is an image showing the downwelling irradiance map that was used for these renderings, with bearing along the horizontal axis and zenith along the vertical:

![image irradiance map](https://github.com/martinvanleeuwen/RenderJay.jl/blob/main/img/power32_ambient0_diffuse10.png)

And here is an overhead view of the Ofenpass scene. You can clearly see that there is a discontinuity along the bearing of the irradiance map -- which does nicely indicate how the map drapes around the visible sky!

![image overhead view](https://github.com/martinvanleeuwen/RenderJay.jl/blob/main/img/Ofenpass_overhead_view_v2.png)



# Example code

Below is an example that renders a set of four Cornell boxes with different colours floating just above a flat surface. If you are on a laptop or you have few logical cores, please mind the addprocs() line and set the number to something comfortable, e.g., the number of logical cores that are available.

```
using Distributed
addprocs(3)
@everywhere using RenderJay
using CSV, LightXML, DataFrames, DelimitedFiles, LinearAlgebra, ProgressMeter, SharedArrays, ImageMagick, Images

pth = pathof(RenderJay)
src_folder, _ = splitdir(pth)
root_folder = src_folder[1:end-3]
scene_fn = joinpath(root_folder, "test/testplot_little_cornellboxes.xml")
assets, geometries, palettes, geometry_bvhs, scene_bvh, skymap, camera = read_scene(scene_fn);
coords = create_coords(camera);
img = SharedArray{Float64, 3}(camera.nBands, camera.xResolution, camera.yResolution);

@time @sync @distributed for coord in coords[1:nprocs()]
    img[:, coord.x, coord.y] = render_pixel(coord, assets, geometries, palettes, geometry_bvhs, scene_bvh, skymap, camera)
end

npixels = length(coords)
blocksize = camera.xResolution
@showprogress for i=1:blocksize:npixels
    e = min(i+blocksize, npixels)
    @sync @distributed for coord in coords[i:e]
        img[:, coord.x, coord.y] = render_pixel(coord, assets, geometries, palettes, geometry_bvhs, scene_bvh, skymap, camera)
    end
end
```

You can save the rendering to disk as follows:

```
using ImageMagick, Images
img3 = reshape(hcat(img...), camera.nBands, camera.xResolution, camera.yResolution)
mn, mx = minimum(img3), maximum(img3)
imscl = img3 ./ (mx/1.0)
imscl .*= 10.0
msk = imscl .> 1.0
imscl[msk] .= 1.0
im = colorview(RGB, imscl)
save("/tmp/test_little_cornell_boxes.tif", im')
a = 1.5
b = 10.0
_, imwb = weibull(imscl, a, b)
im2 = colorview(RGB, imwb)
save("/tmp/test_little_cornell_boxes_wb.tif", im2')
```

Et voila! You should be seeing something like this.

![image Cornell boxes](https://github.com/martinvanleeuwen/RenderJay.jl/blob/main/img/test_little_cornell_boxes_wb_crop.png)

