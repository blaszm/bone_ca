# main program
const n_timesteps = 101
include("_include.jl")

function main()
    if(nprocs() > 1)
        @inbounds for worker in workers()
            rmprocs(worker)
        end
    end

    if(nprocs() == 1)
        addprocs(40)
    end
    @everywhere include("_include.jl")

    # call main with macro- and micromodel, which should be used, and amplitude of applied displacement
    solve_bone_macro("cylinder_air", "rve30pc", 2*1e-6); 
    
    @inbounds for worker in workers()
        rmprocs(worker)
    end
end

#--- simulation parameters
function initmaterialparameters(Δₜ::Float64)
    # cortical bone
    E = 22.e9 # Pa
    ν = 0.32 # -
    ϵ₁ = 8.85e-12 # F/m
    ϵ₃ = 8.85e-12 # F/m
    μᶜ = (1.257e-6)^(-1) # m/H, Inverse!
    eₚ14 = 3.e-3 # As/m^2, d*C 

    mp_b = CorticalBone(E, ν, ϵ₁, ϵ₃, μᶜ, eₚ14)

    # bone marrow
    E = 2.e9 # Pa
    ν = 0.3 # - 
    ϵ₁ = 8.85e-12 # F/m
    ϵ₃ = 8.85e-12 # F/m
    μᶜ = (1.257e-6)^(-1) # m/H, Inverse!
    κ₁ = 1.e4 # 5.e-3 # S/m
    μᵥ = 0.5*1.e-9*Δₜ # s/Pa
    
    mp_m = BoneMarrow(E, ν, ϵ₁, ϵ₃, μᶜ, κ₁, μᵥ)

    # surrounding medium air
    ϵ₁ = 8.85e-12 # F/m
    ϵ₃ = 8.85e-12 # F/m
    μᶜ = (1.257e-6)^(-1) # m/H, Inverse!

    mp_s = SurroundingMediumAir(ϵ₁, ϵ₃, μᶜ)

    return mp_b, mp_m, mp_s
end

# main macroscale solver
function solve_bone_macro(modeltype::String, micromodel::String, amplitude::Float64)
    reset_timer!()
    
    # micro mesh parameters
    a = 0.00032
    b = 0.00036
    n = 2
    
    # preperation export
    if(isdir("Results") == false)
        mkdir("Results")
    end

    sp = SimulationParameters(0.1, 1e-2, 1e-8, 1.0) # ρ_∞, Δₜ, NEWTON_TOL, γ

    # unpack time integration parameters
    α_f = sp.α_f
    αₘ = sp.αₘ
    γₐ = sp.γₐ
    Δₜ = sp.Δₜ
    NEWTON_TOL = sp.NEWTON_TOL*10

    # macromesh - grid, dofhandler, boundary condition
    grid = load_model(modeltype)
    celltype_macro = getcelltype(grid)
    n_dim = Ferrite.getdim(grid.cells[1]) # spacial dimension
    n_npe = Ferrite.nnodes(celltype_macro) # number of nodes per element
    n_fpe = Ferrite.nfaces(celltype_macro) # number of faces per element

    cellvalues_u, _, _ = create_values(celltype_macro)
    dh = create_dofhandler(grid, celltype_macro)
    dbc = create_bc_macro(dh, amplitude)

    # pre-allocate solution vectors, etc.
    n_dofs = ndofs(dh)  # total number of dofs
    n_dofs_n = Int(ndofs_per_cell(dh)/n_npe) # DoF per node
    n_el = getncells(dh.grid) # total number of elements
    n_qp = getnquadpoints(cellvalues_u) # number of QP per element
    
    # current step t_n / start (t=0)
    dn = zeros(n_dofs) # displacement solution vector
    dpn = zeros(n_dofs) # velocity solution vector / time derivative
    vn = zeros(n_dofs) 
    vpn = zeros(n_dofs)

    # new step t_n+1
    dn1 = zeros(n_dofs)
    dpn1 = zeros(n_dofs)
    vn1 = zeros(n_dofs)
    vpn1 = zeros(n_dofs)
    
    # helper variables
    u_r = zeros(n_dofs) # residual displacements
    v_r = zeros(n_dofs) # residual velocities
    a_r = zeros(n_dofs) # residual accelarations
    
    # 3: vpn .= 0 # f(dn, dpn, F0) # calculation necessary, if e.g. F =/= 0 (Eq. 4)
    # 4: vn .= dpn # only for first step - dp/dpn can be =/= 0
    
    Δd = zeros(n_dofs) # solution increment
    r = zeros(n_dofs)  # residual
    S = create_sparsity_pattern(dh); # tangent stiffness matrix
    
    norm_r = 1.0
    norm_r0 = 1.0
    
    # creater helper variables for visualization of different sets
    material_number = zeros(getncells(grid))
    el_mat1 = 0 # number of elements spongybone 
    el_mat2 = 0 # number of elements air
    @inbounds for cell in CellIterator(dh)
        if cellid(cell) in getcellset(grid, "spongybone")
            material_number[cellid(cell)] = 1
            el_mat1 += 1
        elseif cellid(cell) in getcellset(grid, "air")
            material_number[cellid(cell)] = 2
            el_mat2 += 1
        end
    end
    
    cbnumber = zeros(getncells(grid))
    cb_no = 0
    @inbounds for i in 1:length(CellIterator(dh))
        if i in getcellset(grid, "spongybone")
            cbnumber[i] = cb_no
        elseif (i in getcellset(grid, "air"))
            cb_no += 1
        end
    end
    
    # create Material states
    states = [[MacroMaterialState() for _ in 1:n_qp] for _ in 1:getncells(grid)]
    
    # micro material (once)
    mp_b, mp_m, mp_s = initmaterialparameters(Δₜ)
    
    # create micro problem
    mp = MicroProblem(a, b, n, sp, mp_b, mp_m, micromodel)
    celltype_micro = getcelltype(mp.ctx.grid)

    # preparation HDF-files
    hdfname_i = "epsidata_i.h5"
    hdfname_ii = "epsidata_ii.h5"
    n_cells = length(mp.ctx.grid.cells) # n_cells = n_el_RVE
    preparehdf(hdfname_i, hdfname_ii, n_cells, el_mat1, n_qp)
        
    # macro tangent moduli
    mtm = calculatemacrotangent(mp, mp_b, mp_m)
    
    # total Volume 
    Ω = 0.
    @inbounds for (_, cell) in enumerate(CellIterator(dh))
        reinit!(cellvalues_u, cell) 
        @inbounds for GPi in 1:n_qp 
            Ω += getdetJdV(cellvalues_u, GPi)
        end
    end 
    
    # pvd file init 
    pvd = initpvd_macro()
    
    writeinfo_macro(mp_b, mp_m, mp_s, n_dim, n_npe, n_fpe, n_dofs_n, n_qp, n_el, sp, modeltype, micromodel, celltype_macro, celltype_micro)
    
    writebc(dbc, grid, "bc_macro")
        
    # Newton-Raphson loop
    print("\n Starting Newton iterations:\n")
    
    @inbounds for timestep in 1:n_timesteps
        t = timestep # current timestep
        @inbounds dn1[:] = dn[:]
        
        Ferrite.update!(dbc, t)
        apply!(dn1, dbc)

        apply_zero!(dpn, dbc)
        apply_zero!(vn, dbc)
        apply_zero!(vpn, dbc)
        
        newton_itr = -1
        print("\n Time step @time = $timestep:\n")

        while true; newton_itr += 1
            
            if newton_itr > 10
                error("Reached maximum Newton iterations, aborting")
                break
            end
        
            @inbounds vn1 = αₘ/(α_f*γₐ*Δₜ) * (dn1-dn) + (γₐ-αₘ)/(γₐ*α_f) * dpn + (α_f-1)/α_f * vn # (Eq. 18)
            @inbounds vpn1 =  αₘ/(α_f*γₐ^2*Δₜ^2)*(dn1-dn) - vn/(α_f*γₐ*Δₜ) + (γₐ-1)/γₐ * vpn + (γₐ-αₘ)/(α_f*γₐ^2*Δₜ) * dpn # (Eq. 19)

            @inbounds u_r = α_f*dn1+(1-α_f)*dn # (Eq. 9) - dn+α_f
            @inbounds v_r = α_f*vn1+(1-α_f)*vn # (Eq. 10) - vn+α_f
            @inbounds a_r = αₘ*vpn1+(1-αₘ)*vpn # (Eq. 8) - vpn+αₘ
            #@inbounds Fnα_f = α_f*tr_n1 + (1-α_f)*tr_n # (Eq. 11) - necessary if external forces are applied
            
            # assembly and solve
            S, r = doassemble_macro(S, grid, dh, mp_b, mp_m, mp_s, mp, mtm, u_r, v_r, a_r, states, t, cbnumber);

            apply_zero!(S, r, dbc) 
            
            norm_r = norm(r)
            
            if(newton_itr==0 && norm_r > 0.0)
                norm_r0 = norm_r
            end
            
            print("Iteration: $newton_itr \trel. residual: $(@sprintf("%.4e", norm_r/norm_r0))\n")

            if (norm_r/norm_r0) < NEWTON_TOL
                break
            elseif norm_r < NEWTON_TOL
                print("Absolute norm chosen: $(@sprintf("%.4e", norm_r))\n")
                break
            end  
                        
            Δd,_,_,_,_ = KrylovMethods.bicgstb(S, r; maxIter = 1000)
            
            apply_zero!(Δd, dbc)
            @inbounds dn1 .+= Δd

        end # of loop while NR-Iteration
                 
        # compute new accelarations
        @inbounds dpn1 = (dn1 - dn)/(γₐ*Δₜ)+(γₐ-1)/γₐ * dpn # (Eq. 14)
    
        # update t_n+1 -> t_n
        @inbounds dn[:] = dn1[:]
        @inbounds dpn[:] = dpn1[:]
        @inbounds vn[:] = vn1[:]
        @inbounds vpn[:] = vpn1[:] 
        
        # export 
        writeoutput(dh, states, grid, t, dn, vn, vpn, pvd; material_number=material_number)

        updatehdf(hdfname_i, hdfname_ii)
        
    end # of loop timestep  
    
    savepvd(pvd)
    
    print_timer(title = "Sim. $n_timesteps Timesteps", linechars = :ascii)
end;

main();