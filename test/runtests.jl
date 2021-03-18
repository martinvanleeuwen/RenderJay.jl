
using Test

using RenderJay

this_fn = @__FILE__
this_folder, _ = splitdir(this_fn)
scene_fn = joinpath(this_folder, "testplot_little_cornellboxes.xml")

println("Reading data...")
assets, geometries, palettes, geometry_bvhs, scene_bvh, skymap, camera = read_scene(scene_fn)
coords = create_coords(camera)

println("Testing the render_pixel() function...")
radiance = render_pixel(coords[1], assets, geometries, palettes, geometry_bvhs, scene_bvh, skymap, camera)

@test all(radiance .> 0.0)
