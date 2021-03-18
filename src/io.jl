##################################################################
#
#  I/O functionality...
#
##################################################################

# to specify the type of the "geometries" array
allequal(x) = all(y->y==x[1],x)


function read_shaders(shaderXML::Array{XMLElement,1})
    shaders = Function[]
    for el in shaderXML
        n = find_element(el, "function")
        shaderName = content(n)
        s = Symbol(shaderName)
        shader = getfield(Main, s)
        args = child_elements(find_element(el, "args"))
        d = Dict()
        for arg in args
            if occursin(",", content(arg)) # array of values...
                a = parse.(Float64, split(content(arg), ","))
                varName = Symbol(name(arg))
                push!(d, varName => a)
            else # a scalar...
                b = convert(Float64, content(arg))
                varName = Symbol(name(arg))
                push!(d, varName => b)
            end
        end
        # assert that number of arguments equals number of specified spectral bands (i.e., for camera)...
        # not sure how to implement such an assertion...
        f = (Ix, Iy, Iz, vx, vy, vz, nx, ny, nz) -> shader(Ix, Iy, Iz, vx, vy, vz, nx, ny, nz; d...)
        shaders = push!(shaders, f)
        #ds = push!(ds, d)
    end
    return shaders #, ds
end


function read_items(xroot::XMLElement ; add_one=false, load_bvh=true, base_path::AbstractString="")
    if load_bvh
	assets, geometries, palettes, geometry_bvhs = read_items_load_bvh(xroot ; add_one=add_one, base_path=base_path)
        return assets, geometries, palettes, geometry_bvhs
    else
        assets, geometries, palettes = read_items_no_load_bvh(xroot ; add_one=add_one, base_path=base_path)
        return assets, geometries, palettes
    end
end


function read_items_no_load_bvh(xroot::XMLElement ; add_one=false, base_path::AbstractString="")
    oid = 0
    pid = 0
    assets = Asset{Float64, Int64}[]
    geom_args = []
    geometries = Array{Geometry{Float64, Int64, Array{Float64,1}, Array{Int64,1}},1}(undef, 0) # could implement a supertype Geometry?
    palettes = Array{Vector{Function},1}(undef, 0)
    @showprogress for item in xroot["item"]
        # collect geometry...
        geometry_type = attribute(item, "geometry")
        geometry_src = attribute(item, "src")
	if !isabspath(geometry_src)
		geometry_src = joinpath(base_path, geometry_src)
	end
        geom_arg = (geometry_type, geometry_src)
        
        # check if geometry in set; add if not...
        if !(geom_arg in geom_args)
            if geometry_type === "mesh"
                geometry = read_mesh(geometry_src; add_one=add_one, load_bvh=false)
                push!(geometries, geometry)
            elseif geometry_type === "cylinder"
                geometry = read_cylinder(geometry_src; add_one=add_one, load_bvh=false)
                push!(geometries, geometry)
            elseif geometry_type === "cone"
                geometry = read_cone(geometry_src; add_one=add_one, load_bvh=false)
                push!(geometries, geometry)
            elseif geometry_type === "balls"
                geometry = read_ball(geometry_src; add_one=add_one, load_bvh=false)
                push!(geometries, geometry)
            else
                println( string("Cannot read << ", geometry_src, " >> of type << ", geometry_type, " >>, skipping item...") )
                continue
            end
            push!(geom_args, geom_arg)
        end
        
        # find the index...
        oid = -1
        for (j, v) in enumerate(geom_args)
            if geom_arg == v
                oid = j
            end
        end
        
        # collect shaders (notice that palettes do not have to be unique, 
        # unlike geometries which are bigger/consume more memory...
        palette = read_shaders(item["shader"])        
        push!(palettes, palette)
        pid += 1
        
        # collect hooks...
        hooks = item["hook"]
        for hook in hooks
            x = parse(Float64, content(hook["x"][1]))
            y = parse(Float64, content(hook["y"][1]))
            z = parse(Float64, content(hook["z"][1]))
            beta = parse(Float64, content(hook["beta"][1]))
            gamma = parse(Float64, content(hook["gamma"][1]))
            asset = Asset(x, y, z, beta, gamma, oid, pid)
            push!(assets, asset)
        end        
    end

    # conversion may or may not help you, but will first
    # require new functions to be defined for render_pixel,
    # specific for mesh, cylinder, ball, disk...

    ## convert array of type Any[] to a precise composite type, if at all possible...
    #if allequal([typeof(geometries[i]) for i=1:length(geometries)])
    #    geometries = convert.(typeof(geometries[1]), geometries)
    #end
    return assets, geometries, palettes
end


function read_items_load_bvh(xroot::XMLElement ; add_one=false, base_path::AbstractString="")
    oid = 0
    pid = 0
    assets = Asset{Float64, Int64}[]
    geom_args = []
    geometries = Array{Geometry{Float64, Int64, Array{Float64,1}, Array{Int64,1}},1}(undef, 0) # could implement a supertype Geometry?
    geometry_bvhs = Array{Bvh{Float64, Int64, Array{Float64,1}, Array{Int64,1}},1}(undef, 0)
    palettes = Array{Vector{Function},1}(undef, 0)
    @showprogress for item in xroot["item"]
        # collect geometry...
        geometry_type = attribute(item, "geometry")
        geometry_src = attribute(item, "src")
	if !isabspath(geometry_src)
		geometry_src = joinpath(base_path, geometry_src)
	end
        geom_arg = (geometry_type, geometry_src)
        
        # check if geometry in set; add if not...
        if !(geom_arg in geom_args)
            if geometry_type === "mesh"
                geometry, bvh = read_mesh(geometry_src; add_one=add_one, load_bvh=true)
                push!(geometries, geometry)
                push!(geometry_bvhs, bvh)
            elseif geometry_type === "cylinder"
                geometry, bvh = read_cylinder(geometry_src; add_one=add_one, load_bvh=true)
                push!(geometries, geometry)
                push!(geometry_bvhs, bvh)
            elseif geometry_type === "cone"
                geometry, bvh = read_cone(geometry_src; add_one=add_one, load_bvh=true)
                push!(geometries, geometry)
                push!(geometry_bvhs, bvh)
            elseif geometry_type === "balls"
                geometry, bvh = read_ball(geometry_src; add_one=add_one, load_bvh=true)
                push!(geometries, geometry)
                push!(geometry_bvhs, bvh)
            else
                println( string("Cannot read << ", geometry_src, " >> of type << ", geometry_type, " >>, skipping item...") )
                continue
            end
            push!(geom_args, geom_arg)
        end
        
        # find the index...
        oid = -1
        for (j, v) in enumerate(geom_args)
            if geom_arg == v
                oid = j
            end
        end
        
        # collect shaders (notice that palettes do not have to be unique, 
        # unlike geometries which are bigger/consume more memory...
        palette = read_shaders(item["shader"])        
        push!(palettes, palette)
        pid += 1
        
        # collect hooks...
        hooks = item["hook"]
        for hook in hooks
            x = parse(Float64, content(hook["x"][1]))
            y = parse(Float64, content(hook["y"][1]))
            z = parse(Float64, content(hook["z"][1]))
            beta = parse(Float64, content(hook["beta"][1]))
            gamma = parse(Float64, content(hook["gamma"][1]))
            asset = Asset(x, y, z, beta, gamma, oid, pid)
            push!(assets, asset)
        end        
    end

    # conversion may or may not help you, but will first
    # require new functions to be defined for render_pixel,
    # specific for mesh, cylinder, ball, disk...

    # convert array of type Any[] to a precise composite type, if at all possible...
    #if allequal([typeof(geometries[i]) for i=1:length(geometries)])
    #    geometries = convert.(typeof(geometries[1]), geometries)
    #end
    return assets, geometries, palettes, geometry_bvhs
end


function read_camera(xroot::XMLElement)
    u = xroot["camera"][1]
    eyeX = parse(Float64, content(find_element(u, "eyeX")))
    eyeY = parse(Float64, content(find_element(u, "eyeY")))
    eyeZ = parse(Float64, content(find_element(u, "eyeZ")))

    lookX = parse(Float64, content(find_element(u, "lookX")))
    lookY = parse(Float64, content(find_element(u, "lookY")))
    lookZ = parse(Float64, content(find_element(u, "lookZ")))

    eyePoint = [eyeX, eyeY, eyeZ]
    lookAtPoint = [lookX, lookY, lookZ]

    fov = parse(Float64, content(find_element(u, "fov")))
    xRes = parse(Int, content(find_element(u, "xResolution")))
    yRes = parse(Int, content(find_element(u, "yResolution")))
    nBands = parse(Int, content(find_element(u, "nBands")))
    rppx = parse(Int, content(find_element(u, "rppx")))
    nBounces = parse(Int, content(find_element(u, "nBounces")))

    camera = Camera(eyePoint, lookAtPoint, fov, xRes, yRes, nBands, rppx, nBounces)
    return camera
end


function read_skymap(xroot::XMLElement ; base_path::AbstractString="")
    sky_src = attribute(xroot["sky"][1], "src")
    if !isabspath(sky_src)
        sky_src = joinpath(base_path, sky_src)
    end
    skymap = readdlm(sky_src, ',', Float64)
    return skymap
end


function read_scene(scenefn::AbstractString ; add_one=false)
    # scenefn is the full path (and filename) of the scene (XML) file
    abs_scenefn = abspath(scenefn)
    base_path, _ = splitdir(abs_scenefn)
    f = read(abs_scenefn, String)
    xml = parse_file(abs_scenefn)
    xroot = root(xml)
    assets, geometries, palettes, geometry_bvhs = read_items(xroot ; add_one=add_one, load_bvh=true, base_path=base_path)
    camera = read_camera(xroot)
    skymap = read_skymap(xroot ; base_path=base_path)
    abs_scene_bvhfn = string(abs_scenefn[1:end-3], "csv")
    if !isfile(abs_scene_bvhfn)
        scene_bvh = make_bvh(assets, geometry_bvhs)
        write_bvh(scene_bvh, abs_scene_bvhfn)
    end
    scene_bvh = read_bvh(abs_scene_bvhfn)
    return assets, geometries, palettes, geometry_bvhs, scene_bvh, skymap, camera
end


function read_cylinder_data(cylinder_src; add_one=false)
    # loading the axes and diameters...
    mfn = joinpath(cylinder_src, "m.csv")
    lfn = joinpath(cylinder_src, "l.csv")
    dfn = joinpath(cylinder_src, "d.csv")
    mtlfn = joinpath(cylinder_src, "mtl.csv")
    m::Array{Float64, 2} = readdlm(mfn, ',', Float64)
    l::Array{Int, 2} = readdlm(lfn, ',', Int)
    d::Array{Float64, 1} = vec(readdlm(dfn, ',', Float64))
    mtl::Array{Int, 1} = vec(readdlm(mtlfn, ',', Int))
    if add_one
        t += 1
    end
    return m, l, d, mtl
end


function read_cylinder(cylinder_src::AbstractString; add_one=false, load_bvh=true)
    m, l, d, mtl = read_cylinder_data(cylinder_src; add_one=add_one)
    cylinder = Cylinder(m, l, d, mtl)
    if load_bvh === true
        bvhfn = joinpath(cylinder_src, "bvh.csv")
        if !isfile(bvhfn)
            bvh = make_bvh(cylinder)
            write_bvh(bvh, bvhfn)
        end        
	bvh = read_bvh(bvhfn)
        return cylinder, bvh
    else
        return cylinder
    end
end


function read_cone_data(cone_src::AbstractString; add_one=false)
    # loading the axes and diameters...
    mfn = joinpath(cone_src, "m.csv")
    lfn = joinpath(cone_src, "l.csv")
    d1fn = joinpath(cone_src, "d1.csv")
    d2fn = joinpath(cone_src, "d2.csv")
    mtlfn = joinpath(cone_src, "mtl.csv")
    m::Array{Float64, 2} = readdlm(mfn, ',', Float64)
    l::Array{Int, 2} = readdlm(lfn, ',', Int)
    d1::Array{Float64, 1} = vec(readdlm(df1n, ',', Float64))
    d2::Array{Float64, 1} = vec(readdlm(df2n, ',', Float64))
    mtl::Array{Int, 1} = vec(readdlm(mtlfn, ',', Int))
    if add_one
        t += 1
    end
    return m, l, d1, d2, mtl
end

function read_cone(cone_src::AbstractString; add_one=false, load_bvh=true)
    m, l, d, mtl = read_cone_data(cone_src; add_one=add_one)
    cone = Cone(m, l, d, mtl)
    if load_bvh === true
        bvhfn = joinpath(cone_src, "bvh.csv")
        if !isfile(bvhfn)
            bvh = make_bvh(cone)
            write_bvh(bvh, bvhfn)
        end
	bvh = read_bvh(bvhfn)
        return cone, bvh
    else
        return cone
    end
end

function read_disk_data(disk_scr::AbstractString; add_one=false)
    # loading the centers, normal vectors, and diameters...
    pfn = joinpath(disk_scr, "/p.csv")
    nfn = joinpath(disk_scr, "/n.csv")
    dfn = joinpath(disk_scr, "/d.csv")
    mtlfn = joinpath(disk_scr, "/mtl.csv")
    p::Array{Float64, 2} = readdlm(pfn, ',', Float64)
    n::Array{Int, 2} = readdlm(nfn, ',', Int)
    d::Array{Float64, 1} = vec(readdlm(dfn, ',', Float64))
    mtl::Array{Int, 1} = vec(readdlm(mtlfn, ',', Int))
    if add_one
        t += 1
    end
    return p, n, d, mtl
end

function read_disk(disk_scr::AbstractString; add_one=false, load_bvh=true)
    p, n, d, mtl = read_disk_data(disk_scr; add_one=add_one)
    disk = Disk(p, n, d, mtl)
    if load_bvh === true
        bvhfn = joinpath(disk_scr, "/bvh.csv")
        if !isfile(bvhfn)
            bvh = make_bvh(disk)
            write_bvh(bvh, bvhfn)
        end
	bvh = read_bvh(bvhfn)
        return disk, bvh
    else
        return disk
    end
end

function read_ball_data(ball_scr::AbstractString; add_one=false)
    # loading the centers, normal vectors, and diameters...
    pfn = joinpath(ball_scr, "/p.csv")
    dfn = joinpath(ball_src, "/d.csv")
    mtlfn = string(ball_scr, "/mtl.csv")
    p::Array{Float64, 2} = readdlm(pfn, ',', Float64)
    n::Array{Int, 2} = readdlm(nfn, ',', Int)
    d::Array{Float64, 1} = vec(readdlm(dfn, ',', Float64))
    mtl::Array{Int, 1} = vec(readdlm(mtlfn, ',', Int))
    if add_one
        t += 1
    end
    return p, d, mtl
end

function read_ball(ball_scr::AbstractString; add_one=false, load_bvh=true)
    p, d, mtl = read_ball_data(ball_scr; add_one=add_one)
    ball = Ball(p, d, mtl)
    if load_bvh === true
        bvhfn = joinpath(ball_scr, "/bvh.csv")
        if !isfile(bvhfn)
            bvh = make_bvh(ball)
            write_bvh(bvh, bvhfn)
        end
	bvh = read_bvh(bvhfn)
        return ball, bvh
    else
        return ball
    end
end

function read_mesh_data(mesh_src; add_one=false)
    # loading the triangles and creating the mesh...
    vfn = joinpath(mesh_src, "v.csv")
    tfn = joinpath(mesh_src, "t.csv")
    mtlfn = joinpath(mesh_src, "mtl.csv")
    v::Array{Float64, 2} = readdlm(vfn, ',', Float64)
    t::Array{Int, 2} = readdlm(tfn, ',', Int)
    if add_one
        t += 1
    end
    mtl::Array{Int, 1} = vec(readdlm(mtlfn, ',', Int))
    return v, t, mtl
end


function read_mesh(mesh_src::AbstractString; add_one=false, load_bvh=true)
    # mesh_src is the full path of the mesh folder, where the v, t and mtl files are stored
    v, t, mtl = read_mesh_data(mesh_src; add_one=add_one)
    mesh = Mesh(v, t, mtl)
    if load_bvh === true
        bvhfn = joinpath(mesh_src, "bvh.csv")
        if !isfile(bvhfn)
            bvh = make_bvh(mesh)
            write_bvh(bvh, bvhfn)
        end
	bvh = read_bvh(bvhfn)
        return mesh, bvh
    else
        return mesh
    end
end


function write_bvh(bvh::Bvh, bvhfn::AbstractString)
    # bvhfn is the full path (and filename) of the BVH
    df = DataFrame(xmin=bvh.xmin, ymin=bvh.ymin, zmin=bvh.zmin, xmax=bvh.xmax, ymax=bvh.ymax, zmax=bvh.zmax, bvh2bb=bvh.bvh2bb, left_child=bvh.left_child, right_child=bvh.right_child)
    CSV.write(bvhfn, df, delim=" ")
end


function read_bvh(bvhfn::AbstractString)
    # bvhfn is the full path (and filename) of the BVH
    df = CSV.read(bvhfn, DataFrame) #delim=" ")
    xmin::Array{Float64}, ymin::Array{Float64}, zmin::Array{Float64}, xmax::Array{Float64}, ymax::Array{Float64}, zmax::Array{Float64} = convert(Array{Float64}, df.xmin), convert(Array{Float64}, df.ymin), convert(Array{Float64}, df.zmin), convert(Array{Float64}, df.xmax), convert(Array{Float64}, df.ymax), convert(Array{Float64}, df.zmax)
    bvh2bb::Array{Int}, right_child::Array{Int}, left_child::Array{Int} = convert(Array{Int}, df.bvh2bb), convert(Array{Int}, df.right_child), convert(Array{Int}, df.left_child)
    bvh = Bvh{Float64, Int, Array{Float64,1}, Array{Int,1}}(xmin, ymin, zmin, xmax, ymax, zmax, bvh2bb, left_child, right_child)
    return Bvh(xmin, ymin, zmin, xmax, ymax, zmax, bvh2bb, left_child, right_child)
end


