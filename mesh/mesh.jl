# function to set up cell-/facevalues- and facevalues (choose interpolation rules, i.e. order of shape functions
# for geometry and potential fields and quadrature points for numerical integration
function create_values(celltype)

    # QuadratureRule{dim,Shape}(Order) : order 2 is standard integration scheme with +/- 1/sqrt3 quadrature points
    # geometric interpolation: Lagrange{}(dim,Shape,Order) : order 1 means linear shape functions 

    if(celltype == Hexahedron)
        qr = QuadratureRule{3,RefCube}(2)
        interpolation_geom = Lagrange{3,RefCube,1}()
        cellvalues_u = CellVectorValues(qr, Lagrange{3,RefCube,1}(), interpolation_geom)
        cellvalues_phi = CellScalarValues(qr, Lagrange{3,RefCube,1}(), interpolation_geom)
        cellvalues_A = CellVectorValues(qr, Lagrange{3,RefCube,1}(), interpolation_geom)
        return cellvalues_u, cellvalues_phi, cellvalues_A
    elseif(celltype == Tetrahedron)
        qr = QuadratureRule{3,RefTetrahedron}(2)
        interpolation_geom = Lagrange{3,RefTetrahedron,1}()
        cellvalues_u = CellVectorValues(qr, Lagrange{3,RefTetrahedron,1}(), interpolation_geom)
        cellvalues_phi = CellScalarValues(qr, Lagrange{3,RefTetrahedron,1}(), interpolation_geom)
        cellvalues_A = CellVectorValues(qr, Lagrange{3,RefTetrahedron,1}(), interpolation_geom)
        return cellvalues_u, cellvalues_phi, cellvalues_A
    end

end;

# function to set up the degrees of freedom 
function create_dofhandler(grid::Grid, celltype)
    dh = DofHandler(grid)
    # u: mechanical displacement (3 DoF)
    # phi: electric scalar potential (1 DoF)
    # A: magnetic vector potential (3 DoF)
    if(celltype == Hexahedron)
        push!(dh, :u, 3, Lagrange{3,RefCube,1}())
        push!(dh, :phi, 1, Lagrange{3,RefCube,1}())
        push!(dh, :A, 3, Lagrange{3,RefCube,1}())
    elseif(celltype == Tetrahedron)
        push!(dh, :u, 3, Lagrange{3,RefTetrahedron,1}()) 
        push!(dh, :phi, 1, Lagrange{3,RefTetrahedron,1}()) 
        push!(dh, :A, 3, Lagrange{3,RefTetrahedron,1}()) 
    end
    close!(dh)
    return dh
end;

# kann raus sobald PBC über Ferrite direkt statt über CoheStr implementiert sind
function create_dofhandler_helper(grid::Grid, celltype)
    dh = DofHandler(grid)
    if(celltype == Hexahedron)
        push!(dh, :T, 1, Lagrange{3,RefCube,1}())
    elseif(celltype == Tetrahedron) 
        push!(dh, :T, 1, Lagrange{3,RefTetrahedron,1}())
    end
    close!(dh)
    return dh
end; 

# microscale boundary conditions: fix corner nodes completely
function create_bc(dh::DofHandler) 
    dbc = ConstraintHandler(dh)
    add!(dbc, Dirichlet(:u, getnodeset(dh.grid, "corners"), (x,t) -> zero(Vec{3}), [1,2,3]))
    add!(dbc, Dirichlet(:phi, getnodeset(dh.grid, "corners"), (x,t) -> zero(Vec{1}), [1])) 
    add!(dbc, Dirichlet(:A, getnodeset(dh.grid, "corners"), (x,t) -> zero(Vec{3}), [1,2,3])) 
    close!(dbc)
    t = 0.0 # t = time
    Ferrite.update!(dbc, t)
    return dbc
end;