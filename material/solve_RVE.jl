# microscale solver
function solve_RVE(mq::MacroQuantities, mp::MicroProblem, mp_b::CorticalBone, mp_m::BoneMarrow, globgpnumber::Int; overwrite::Bool=true, output::Bool=false)

    # unpack parameters
    ε_macro = mq.ε_macro
    E_macro = mq.E_macro
    B_macro = mq.B_macro
    #t_macro = mq.t

    a = mp.a
    b = mp.b
    n = mp.n
    ctx = mp.ctx
    bd = mp.bd
    S_red = mp.S_red
    
    α_f = mp.sp.α_f
    αₘ = mp.sp.αₘ
    γₐ = mp.sp.γₐ
    Δₜ = mp.sp.Δₜ  
    NEWTON_TOL = mp.sp.NEWTON_TOL
    nodedoflist = mp.nodedoflist
    
    # HDF file names
    hdfname_i = "epsidata_i.h5"
    hdfname_ii = "epsidata_ii.h5"
    dataname = "epsi_" * string(globgpnumber)

    n_timesteps_micro = 1 
    
    # grid, dofhandler, boundary conditions
    celltype = getcelltype(ctx.grid)
    dh = create_dofhandler(ctx.grid, celltype)
    cellvalues_u, cellvalues_phi, cellvalues_A = create_values(celltype)
    dbc = create_bc(dh)
    
    # pre-allocate solution vectors, etc.
    n_dim = Ferrite.getdim(ctx.grid.cells[1]) # number of dimensions
    n_dofs = ndofs(dh)  # total number of dofs
    n_dofs_n = Int(ndofs_per_cell(dh)/Ferrite.nnodes(celltype)) # DoF per node
    n_el = getncells(dh.grid) # total number of elements
    n_qp = getnquadpoints(cellvalues_u) # number of quadrature points per element
    
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
    
    # 3
    # vpn .= 0 # f(dn, dpn, F0) # calculation necessary, if e.g. F =/= 0 (Eq. 4)
    # 4
    # vn .= dpn # only in first step - dp/dpn can be =/= 0
    
    r = zeros(n_dofs) # residual

    apply!(dn, dbc)
    
    r_red = ones(Int(nBCDofs(ctx, bd)*n_dofs_n))
    norm_r = 1.0
    norm_r0 = 1.0
    
    # creater helper variables for visualization of different sets
    material_number = zeros(getncells(ctx.grid))
    el_mat1 = 0 # number of elements Material_1 (cortical bone)
    el_mat2 = 0 # number of elements Material_2 (bone marrow)
    @inbounds for cell in CellIterator(dh)
        if cellid(cell) in getcellset(ctx.grid, "Material_1")
            material_number[cellid(cell)] = 1
            el_mat1 += 1
        elseif cellid(cell) in getcellset(ctx.grid, "Material_2")
            material_number[cellid(cell)] = 2
            el_mat2 += 1
        end
    end
    
    # create material states
    states = [[MaterialState() for _ in 1:n_qp] for _ in 1:getncells(ctx.grid)] 
    
    # read εⁱ from HDF file and pass to states
    εⁱ_hdf = readhdf(hdfname_i, dataname)
    @inbounds for (elmt, _) in enumerate(CellIterator(dh)) 
        state = states[elmt]
        @inbounds for QPi in 1:n_qp
            j = n_qp*(elmt-1) + QPi
            state[QPi].εⁱ[:] = εⁱ_hdf[:,j]
        end
    end
    
    # total Volume 
    Ω = 0.0
    @inbounds for (_, cell) in enumerate(CellIterator(dh))
        reinit!(cellvalues_u, cell) 
        @inbounds for QPi in 1:n_qp 
            Ω += getdetJdV(cellvalues_u, QPi) 
        end
    end  
    
    # output variables
    σ_aver = zeros(6)
    D_aver = zeros(3)
    Dp_aver = zeros(3)
    H_aver = zeros(3)
    J_aver = zeros(3) 
    
    # pvd file init
    if(output)
        pvd = initpvd(globgpnumber)
    end
    
    # solver
    if(output)
        print("\n Starting Newton iterations:\n")
    end

    for timestep in 1:n_timesteps_micro # time loop
        
        t = timestep # current timestep
        newton_itr = -1
        if(output)
            print("\n Time step @time = $timestep:\n")
        end
        
        @inbounds dn1[:] = dn[:]
        
        while true; newton_itr += 1
            
            if newton_itr > 10
                error("Reached maximum Newton iterations, aborting")
                break
            end
        
            @inbounds vn1 = αₘ/(α_f*γₐ*Δₜ) * (dn1-dn) + (γₐ-αₘ)/(γₐ*α_f) * dpn + (α_f-1)/α_f * vn # (Eq. 18)
            @inbounds vpn1 =  αₘ/(α_f*γₐ^2*Δₜ^2)*(dn1-dn) - vn/(α_f*γₐ*Δₜ) + (γₐ-1)/γₐ * vpn + (γₐ-αₘ)/(α_f*γₐ^2*Δₜ) * dpn # (Eq. 19)

            @inbounds u_r = α_f*dn1+(1-α_f)*dn # (Gl.9) - dn+α_f
            @inbounds v_r = α_f*vn1+(1-α_f)*vn # (Gl.10) - vn+α_f
            @inbounds a_r = αₘ*vpn1+(1-αₘ)*vpn # (Gl.8) - vpn+αₘ
            # use (Eq. 11) if external forces are applied
            
            # assembly and solve
            r = doassemble_r(cellvalues_u, cellvalues_phi, cellvalues_A, ctx.grid, dh, mp_b, mp_m, mp.sp, u_r, v_r, a_r, states, mq); 

            apply_zero!(r, dbc)
            
            Δd_red = zeros(Int(nBCDofs(ctx, bd)*n_dofs_n))

            r_red = applyPBC_r(ctx, bd, r, nodedoflist, n_dofs_n)

            norm_r = norm(r_red)
            
            if(newton_itr==0 && norm_r > 0.0)
                norm_r0 = norm_r
            end
            
            if(output)
                print("Iteration: $newton_itr \trel. residual: $(@sprintf("%.4e", norm_r/norm_r0))\n")
            end

            if (norm_r/norm_r0) < NEWTON_TOL
                break
            elseif (norm_r < NEWTON_TOL && newton_itr > 0)
                if(output)
                    print("Absolute norm chosen: $(@sprintf("%.4e", norm_r))\n")
                end
                break
            end  

            Δd_red,_,_,_,_ = KrylovMethods.bicgstb(S_red, r_red; maxIter = 100)

            Δd = undoPBC(ctx, bd, n_dofs, n_dofs_n, Δd_red, nodedoflist)
            apply_zero!(Δd, dbc)
            @inbounds dn1 .+= Δd

        end # of loop while NR-Iteration
            
        # update all the material states (including update inelastic strain (internal variable)) after equilibrium is reached
        for cell_states in states
            foreach(update_state!, cell_states) 
        end
        
        # compute new accelarations
        @inbounds dpn1 = (dn1 - dn)/(γₐ*Δₜ)+(γₐ-1)/γₐ * dpn # (Eq. 14)
    
        # (Eq. 17) update t_n+1 -> t_n
        @inbounds dn[:] = dn1[:]
        @inbounds dpn[:] = dpn1[:]
        @inbounds vn[:] = vn1[:]
        @inbounds vpn[:] = vpn1[:] 
        
        # volume averages
        σ_aver, D_aver, Dp_aver, H_aver, J_aver = volaver(cellvalues_A, dh, n_qp, Ω, states)
    
        # export
        if(output)
            writeoutput(dh, states, ctx.grid, t, dn, vn, vpn, pvd; ismicro=true, globgpnumber=globgpnumber, material_number=material_number)
        end
            
    end # of loop timestep  
        
    if(output)
        # write info.txt
        writeinfo_micro(mp_b, mp_m, celltype, n_dim, n_dofs_n, n_qp, n_el, Ω, ε_macro, E_macro, B_macro, σ_aver, D_aver, Dp_aver, H_aver, J_aver, el_mat1, el_mat2, globgpnumber, a, b, n)
        # write boundary_conditions.vtu
        writebc(dbc, ctx.grid, "bc_micro")
    	savepvd(pvd)
    end
    
    # collect and save new εⁱ in HDF file
    if(overwrite)
        @inbounds for (elmt, _) in enumerate(CellIterator(dh)) 
            state = states[elmt]
            @inbounds for QPi in 1:n_qp
                j = 8*(elmt-1) + QPi
                εⁱ_hdf[:,j] = state[QPi].εⁱ[:]
            end
        end
        overwritehdf(hdfname_ii, dataname, εⁱ_hdf)
    end
    
    return σ_aver, D_aver, Dp_aver, H_aver, J_aver
end