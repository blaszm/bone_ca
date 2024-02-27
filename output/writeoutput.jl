# function to write output data to vtk files - use kwargs to differentiate between macro (default) and micro output
# important note: for projection of scalar/vector quantities, use project_to_nodes=true, for tensor quantities use project_to_nodes=false - will be removed in new Ferrite release
function writeoutput(dh::DofHandler, states, grid::Grid, t::Int, u::Array{Float64,1}, v::Array{Float64,1}, a::Array{Float64,1}, pvd::WriteVTK.CollectionFile; ismicro::Bool=false, globgpnumber::Int=1, material_number::Array{Float64,1}=zeros(1))
    
    work_dir = pwd()
    cd("Results")

    # Output stresses / strain - x,y,z,xy,yz,xz
    sigma_qp = collect_tensor_quantity(dh, states, :σ)
    eps_qp = collect_tensor_quantity(dh, states, :ε)
       
    # Output other fields    
    D_qp = collect_vector_quantity(dh, states, :D)
    Dp_qp = collect_vector_quantity(dh, states, :Dp)
    E_qp = collect_vector_quantity(dh, states, :E)
    J_qp = collect_vector_quantity(dh, states, :J)
    H_qp = collect_vector_quantity(dh, states, :H)
    B_qp = collect_vector_quantity(dh, states, :B)
    
    # conversion quadratures points -> nodes
    celltype = getcelltype(grid)
    qr = nothing
    projector = nothing
    if(celltype == Hexahedron)
        qr = QuadratureRule{3,RefCube}(2)
        projector = L2Projector(Lagrange{3,RefCube,1}(), grid)
    elseif(celltype == Tetrahedron)
        qr = QuadratureRule{3,RefTetrahedron}(2)
        projector = L2Projector(Lagrange{3,RefTetrahedron,1}(), grid)
    end

    sigma_nodes = project(projector, sigma_qp, qr; project_to_nodes=false)
    eps_nodes = project(projector, eps_qp, qr; project_to_nodes=false)
        
    D_nodes = project(projector, D_qp, qr) 
    Dp_nodes = project(projector, Dp_qp, qr) 
    E_nodes = project(projector, E_qp, qr) 
    J_nodes = project(projector, J_qp, qr) 
    H_nodes = project(projector, H_qp, qr)
    B_nodes = project(projector, B_qp, qr) 
        
    filename = ""
    time = @sprintf("%3.3i", t)
    if(ismicro)
        gpnumber = @sprintf("%5.5i", globgpnumber)
        filename = string("micro_", gpnumber, "_", time)
    else
        filename = string("macro_", time)
    end

    vtk_grid(filename, grid; compress=false) do vtk # creates Paraview output
        vtk_point_data(vtk, dh, u) 
        vtk_point_data(vtk, dh, v, "p") 
        vtk_point_data(vtk, dh, a, "pp")
        
        vtk_point_data(vtk, projector, sigma_nodes, "sigma")
        vtk_point_data(vtk, projector, eps_nodes, "eps")

        # Output inelastic strain (x,y,z,xy,yz,xz) and material number - only for micro problem
        if(ismicro)
            epsi_qp = collect_tensor_quantity(dh, states, :εⁱ)
            epsi_nodes = project(projector, epsi_qp, qr; project_to_nodes=false)
            vtk_point_data(vtk, projector, epsi_nodes, "eps_i")
        end
            
        vtk_point_data(vtk, D_nodes, "D")
        vtk_point_data(vtk, Dp_nodes, "Dp")
        vtk_point_data(vtk, E_nodes, "E")
        vtk_point_data(vtk, J_nodes, "J")
        vtk_point_data(vtk, H_nodes, "H")
        vtk_point_data(vtk, B_nodes, "B")
            
        vtk_cell_data(vtk, material_number, "Mat")

        pvd[t-1] = vtk
    end # of export
    
    cd(work_dir)
    return
end

# visualize boundary conditions
function writebc(dbc, grid::Grid, filename::String)
    work_dir = pwd()
    cd("Results")
    vtk_grid(filename, grid; compress=false) do vtk
        vtk_point_data(vtk, dbc)
    end
    cd(work_dir)
    return
end

# functions for pvd collections
function initpvd(globgpnumber)
    work_dir = pwd()
    cd("Results")
    gpnumber = @sprintf("%5.5i", globgpnumber)
    filename = string("pvd_micro_", gpnumber)
    pvd = paraview_collection(filename)
    cd(work_dir)
    return pvd
end

function initpvd_macro()
    work_dir = pwd()
    cd("Results")
    filename = "pvd_macro"
    pvd = paraview_collection(filename)
    cd(work_dir)
    return pvd
end

function savepvd(pvd)
    work_dir = pwd()
    cd("Results")
    vtk_save(pvd)
    cd(work_dir)
    return
end