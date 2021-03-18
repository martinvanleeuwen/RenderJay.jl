function load_scene_data(scene_full_path::AbstractString)
    labels = Int[]
    X = Float64[]
    Y = Float64[]
    Z = Float64[]
    B = Float64[]
    G = Float64[]
    L = Int[]
    asset_names = String[]
    skymap_name = "nothing"
    assets_fn = "nothing"
    optic_name = "nothing"
    println( string("Opening scene file: ", scene_full_path, " ..."))
    open(scene_full_path) do f
      for line in eachline(f)
        if length(line) > 1
          if occursin("//", line)
            continue
          elseif occursin("#", line)
            label, asset_name = split(line, ",")
            labels = push!(labels, parse(Int, label[2:end]) )
            asset_names = push!(asset_names, asset_name)
          elseif occursin("@", line)
            assets_fn = split(line, ",")[2]
          elseif occursin("*", line)
            skymap_name = split(line, ",")[2]
          elseif occursin("&", line)
            optic_name = split(line, ",")[2]
          else
            x_, y_, z_, beta_, gamma_, label_ = split(line, ",")
            x, y, z = parse(Float64, x_), parse(Float64, y_), parse(Float64, z_)
            beta, gamma = parse(Float64, beta_), parse(Float64, gamma_)
            l = round(Int, parse(Float64, label_))
            X = push!(X, x)
            Y = push!(Y, y)
            Z = push!(Z, z)
            B = push!(B, beta)
            G = push!(G, gamma)
            L = push!(L, l)
          end
        end
      end
    end
    return skymap_name, optic_name, assets_fn, labels, asset_names, L, X, Y, Z, B, G
end

function read_scene(scene_name::AbstractString)
    scene_full_path = string(scenefolder, scene_name, ".jay")
    skymap_name, optic_name, assets_fn, asset_labels, asset_names, asset_idx, x, y, z, beta, gamma = load_scene_data(scene_full_path)

    # load assets...
    xml_path = string(assetfolder, assets_fn, ".xml")
    xdoc = parse_file(xml_path)
    xroot = LightXML.root(xdoc)
    unique_mesh_names, mesh_idx, unique_palette_names, palette_idx = get_assets(xroot, asset_names)
    asset2mesh_idx = mesh_idx[asset_idx]
    asset2palette_idx = palette_idx[asset_idx]
    n = length(x)
    assets = Asset{Float64, Int}[]
    for i=1:n
        push!(assets, Asset(x[i], y[i], z[i], beta[i], gamma[i], asset2mesh_idx[i], asset2palette_idx[i]))
    end
    
    # load mesh data...
    meshes::Array{Mesh{Float64, Int, Array{Float64,1}, Array{Int64,1}},1} = []
    mesh_bvhs::Array{Bvh{Float64, Int, Array{Float64,1}, Array{Int64,1}},1} = []
    for mesh_name in unique_mesh_names
        mesh = read_mesh(mesh_name)
        bvh = read_bvh(mesh_name)
        meshes = push!(meshes, mesh)
        mesh_bvhs = push!(mesh_bvhs, bvh)
    end

    # load scene bvh (or make it if necessary)...
    if !isfile( string(bvhfolder, scene_name, ".bvh") )
        make_bvh(assets, meshes, mesh_bvhs)
    end
    scene_bvh::Bvh{Float64, Int, Array{Float64,1}, Array{Int,1}} = read_bvh(scene_name)

    # load the camera...
    camera = read_camera(optic_name)

    # read the skymap (optional, only if specified in scene description file)...
    skymap_path = string(skyfolder, skymap_name, ".csv")
    skymap = readdlm(skymap_path, ',', Float64)

    # load materials/shaders...
    palettes = Vector{Function}[]
    for palette_name in unique_palette_names
        xmlfn = string(palettefolder, palette_name, ".xml")
        shaders = load_shaders(xmlfn ; assert_nBands=camera.nBands) # shaders, or palette, is of type: Vector{Function}
        palettes = push!(palettes, shaders)
    end

    return assets, meshes, palettes, mesh_bvhs, scene_bvh, skymap, camera
end
