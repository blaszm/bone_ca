# struct which stores all important information for the microscale problem
struct MicroProblem # remove ctx, bd, nodedoflist, add grid instead
    a::Float64 # RVE length
    b::Float64 # RVE length
    n::Int64 # RVE mesh resolution
    sp::SimulationParameters
    ctx::CoherentStructures.GridContext # RVE-mesh and utilities <- kann auch nochmal problematisch werden!
    bd::BoundaryData
    S_red::SparseMatrixCSC{Float64,Int64} # stiffness matrix RVE
    nodedoflist::Array{Int64,2}
end

# constructor for struct MicroProblem - compute only once for all RVEs
function MicroProblem(a::Float64, b::Float64, n, sp::SimulationParameters, mp_b::CorticalBone, mp_m::BoneMarrow, rvemesh::String) 
    #create_rve_msh(a, b, n) # check if file already exists
    
    ctx = nothing
    loc = nothing
    n_npe = nothing

    meshpath = "mesh/"*rvemesh*".msh"

    grid = FerriteGmsh.togrid(meshpath)

    celltype = getcelltype(grid)

    addfaceset!(grid, "left", x -> x[1] ≈ -a-b/2)
    addfaceset!(grid, "right", x -> x[1] ≈ a+b/2)
    addfaceset!(grid, "top", x -> x[2] ≈ -a-b/2)
    addfaceset!(grid, "bottom", x -> x[2] ≈ a+b/2)
    addfaceset!(grid, "front", x -> x[3] ≈ -a-b/2)
    addfaceset!(grid, "back", x -> x[3] ≈ a+b/2)
    addnodeset!(grid, "corners", x -> norm(x[1]) ≈ a+b/2 && norm(x[2]) ≈ a+b/2 && norm(x[3]) ≈ a+b/2)

    # preparation PBC -> kann geändert werden ab hier
    ll = Vec{3}((-a-b/2, -a-b/2,-a-b/2))
    ur = Vec{3}((a+b/2, a+b/2, a+b/2))

    if(celltype == Hexahedron)
        loc = CoherentStructures.Regular3DGridLocator{Ferrite.Hexahedron}((n+1), (n+1), (n+1), ll, ur)
    elseif(celltype == Tetrahedron)
        loc = CoherentStructures.Regular3DGridLocator{Ferrite.Tetrahedron}((n+1), (n+1), (n+1), ll, ur)
    end

    dh2 = create_dofhandler_helper(grid, celltype)

    if(celltype == Hexahedron)
        ctx = CoherentStructures.GridContext{3}(grid, Lagrange{3,RefCube,1}(), Lagrange{3,RefCube,1}(), dh2, QuadratureRule{3,RefCube}(2), loc)
        n_npe = 8 # number of nodes per element
    elseif(celltype == Tetrahedron)
        ctx = CoherentStructures.GridContext{3}(grid, Lagrange{3,RefTetrahedron,1}(), Lagrange{3,RefTetrahedron,1}(), dh2, QuadratureRule{3,RefTetrahedron}(2), loc)
        n_npe = 4 # number of nodes per element
    end

    grid = nothing
    
    predicate = (p1, p2) -> peuclidean(p1, p2, [2*a+b, 2*a+b, 2*a+b]) < 1e-8 
    bd = BoundaryData(ctx, predicate)
    
    dh = create_dofhandler(ctx.grid, celltype)
    cellvalues_u, cellvalues_phi, cellvalues_A = create_values(celltype) 

    dbc = create_bc(dh)
    S = create_sparsity_pattern(dh)

    nodedoflist = nodedofs(dh, ctx.grid)

    doassemble_S!(cellvalues_u, cellvalues_phi, cellvalues_A, S, ctx.grid, dh, mp_b, mp_m, sp);
    
    n_dofs = ndofs(dh)
    n_dofs_n = Int(ndofs_per_cell(dh)/n_npe) # DoF per node
    r = zeros(n_dofs) # helper vector
    apply_zero!(S, r, dbc) 
    S_red = applyPBC_S(ctx, bd, S, nodedoflist, n_dofs_n)
    return MicroProblem(a, b, n, sp, ctx, bd, S_red, nodedoflist)
end