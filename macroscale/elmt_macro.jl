# ue = element displacement solution vector (at nodes)
# ve = element velocity solution vector (at nodes)
# ae = element acceleration solution vector (at nodes)
# ue stores for each node first all displacements u, then electric scalar potential phi, then magnetic vector potential A

# macroscopic element (cell) routine
function elmt_macro!(Se::PseudoBlockArray{Float64,2}, re::PseudoBlockArray{Float64,1}, elmtno::Int, dh::DofHandler, cellvalues_u::CellVectorValues{dim}, cellvalues_phi::CellScalarValues{dim}, cellvalues_A::CellVectorValues{dim}, grid::Grid, mp_b::CorticalBone, mp_m::BoneMarrow, mp::MicroProblem, mtm::MacroTangentModuli, ue::Array{Float64,1}, ve::Array{Float64,1}, ae::Array{Float64,1}, state::Array{MacroMaterialState{Float64,Array{Float64,1}},1}, t::Int, cbnumber::Float64) where {dim} 
    
    n_basefuncs_u = getnbasefunctions(cellvalues_u) # number of shape functions
    n_basefuncs_ϕ = getnbasefunctions(cellvalues_phi)
    n_basefuncs_A = getnbasefunctions(cellvalues_A)
    
    n_qp = getnquadpoints(cellvalues_u) # number of quadrature points
    
    u▄, ϕ▄ ,A▄ = 1,2,3 # block index for coupled problems

    coordinates = [zero(Vec{3}) for _ in 1:length(grid.cells[1].nodes)]
    
    nodeids = grid.cells[elmtno].nodes
    @inbounds for j in eachindex(coordinates)
        coordinates[j] = grid.nodes[nodeids[j]].x
    end 
    
    reinit!(cellvalues_u, coordinates) 
    reinit!(cellvalues_phi, coordinates)
    reinit!(cellvalues_A, coordinates)
    
    # unpack simulation parameters
    ctan = zeros(3)
    ctan[:] = mp.sp.ctan[:]
    γ = mp.sp.γ
    
    # macro tangent moduli
    C = mtm.C
    ϵₜ = mtm.ϵₜ
    μⁱ = mtm.μⁱ
    eₚ = mtm.eₚ
    κ = mtm.κ
    
    # main variables (nodes)
    u = zeros(n_basefuncs_u)
    ϕ = zeros(n_basefuncs_ϕ)
    A = zeros(n_basefuncs_A)
    up = zeros(n_basefuncs_u)
    ϕp = zeros(n_basefuncs_ϕ) 
    Ap = zeros(n_basefuncs_A)
    App = zeros(n_basefuncs_A)

    u[:] = ue[dof_range(dh, :u)]  
    ϕ[:] = ue[dof_range(dh, :phi)] 
    A[:] = ue[dof_range(dh, :A)]
    up[:] = ve[dof_range(dh, :u)] 
    ϕp[:] = ve[dof_range(dh, :phi)] 
    Ap[:] = ve[dof_range(dh, :A)] 
    App[:] = ae[dof_range(dh, :A)] 
    
    @inbounds for QPi in 1:n_qp # loop QPi
        
        # volume part
        dΩ = getdetJdV(cellvalues_u, QPi)   

        # strain + time derivative
        ε = tovoigt(symmetric(function_gradient(cellvalues_u, QPi, u)), offdiagscale = 2, order = [1 4 6; 7 2 5; 9 8 3])
        
        # electric field
        E_t = -1.0*function_gradient(cellvalues_phi, QPi, ϕ) - 1.0*function_value(cellvalues_A, QPi, Ap)
        
        # magnetic flux density
        B_t = function_curl(cellvalues_A, QPi, A)      
        
        # --- start scale transition ---
        # convert Tensors.Vec to Base.Vector
        E, B = zeros(3), zeros(3)
        for i in 1:3
            E[i] = E_t[i]
            B[i] = B_t[i]
        end
        mq = MacroQuantities(ε, E, B, t)
        globgpnumber = Int(getnquadpoints(cellvalues_u)*(elmtno-cbnumber-1) + QPi)
        @timeit "solveRVE" σ, D, Dp, H, J = solve_RVE(mq, mp, mp_b, mp_m, globgpnumber) 
        # --- end scale transition --- 
        
        # save material state
        state[QPi].σ = σ 
        state[QPi].ε = ε
        state[QPi].E = E
        state[QPi].D = D
        state[QPi].Dp = Dp
        state[QPi].B = B
        state[QPi].H = H
        state[QPi].J = J
        
        # calculate stiffness / damping / mass matrices

        # Ke_uu += (BuT * C * Bu) * dΩ
        # Ke_uϕ += (BuT * eₚ' * Bgrad) * dΩ
        # Ce_uA += (BuT * eₚ' * N_A) * dΩ
        @inbounds for i in 1:n_basefuncs_u
            BuT = tovoigt(shape_symmetric_gradient(cellvalues_u, QPi, i), offdiagscale = 2, order = [1 4 6; 7 2 5; 9 8 3])'
            @inbounds for j in 1:n_basefuncs_u
                Bu = tovoigt(shape_symmetric_gradient(cellvalues_u, QPi, j), offdiagscale = 2, order = [1 4 6; 7 2 5; 9 8 3])
                Se[BlockIndex((u▄,u▄), (i,j))] += ctan[1]*((BuT * C * Bu) * dΩ)
            end
            @inbounds for j in 1:n_basefuncs_ϕ
                Bgrad = shape_gradient(cellvalues_phi, QPi, j) 
                Se[BlockIndex((u▄,ϕ▄), (i,j))] += ctan[1]*((BuT * eₚ' * Bgrad) * dΩ)
            end
            @inbounds for j in 1:n_basefuncs_A
                N_A = shape_value(cellvalues_A, QPi, j)
                Se[BlockIndex((u▄,A▄), (i,j))] += ctan[2]*((BuT * eₚ' * N_A) * dΩ)
            end
        end

        # Ke_ϕu += (BgradT ⋅ (eₚ * Bu)) * dΩ
        # Ke_ϕϕ += -1.0 * (BgradT ⋅ (ϵₜ * Bgrad)) * dΩ  
        # Ce_ϕA += -1.0 * (BgradT ⋅ (ϵₜ * N_A)) * dΩ
        @inbounds for i in 1:n_basefuncs_ϕ
            BgradT = shape_gradient(cellvalues_phi, QPi, i)
            @inbounds for j in 1:n_basefuncs_u
                Bu = tovoigt(shape_symmetric_gradient(cellvalues_u, QPi, j), offdiagscale = 2, order = [1 4 6; 7 2 5; 9 8 3])
                Se[BlockIndex((ϕ▄,u▄), (i,j))] += ctan[1]*((BgradT ⋅ (eₚ * Bu)) * dΩ)
            end
            @inbounds for j in 1:n_basefuncs_ϕ
                Bgrad = shape_gradient(cellvalues_phi, QPi, j)
                Se[BlockIndex((ϕ▄,ϕ▄), (i,j))] += ctan[1]*(-1.0 * (BgradT ⋅ (ϵₜ * Bgrad)) * dΩ)
            end
            @inbounds for j in 1:n_basefuncs_A
                N_A = shape_value(cellvalues_A, QPi, j)
                Se[BlockIndex((ϕ▄,A▄), (i,j))] += ctan[2]*(-1.0 * (BgradT ⋅ (ϵₜ * N_A)) * dΩ)
            end
        end

        # Ke_AA += (BcurlT ⋅ (μⁱ * Bcurl) + γ * B_divdyad) * dΩ
        @inbounds for i in 1:n_basefuncs_A
            BcurlT = shape_curl(cellvalues_A, QPi, i)
            BdivT = shape_divergence(cellvalues_A, QPi, i)
            @inbounds for j in 1:n_basefuncs_A
                Bcurl = shape_curl(cellvalues_A, QPi, j)
                Bdiv = shape_divergence(cellvalues_A, QPi, j)
                B_divdyad = Bdiv ⋅ BdivT
                Se[BlockIndex((A▄,A▄), (i,j))] += ctan[1]*((BcurlT ⋅ (μⁱ * Bcurl) + γ * B_divdyad) * dΩ)
            end
        end

        # Ce_Au += -1.0 * (N_AT ⋅ (eₚ * Bu)) * dΩ
        # Ke_Aϕ += (N_AT ⋅ κ * Bgrad) * dΩ
        # Ce_Aϕ += (N_AT ⋅ (ϵₜ * Bgrad)) * dΩ
        # Ce_AA += (N_AT ⋅ κ * N_A) * dΩ
        # Me_AA += (N_AT ⋅ (ϵₜ * N_A)) * dΩ
        @inbounds for i in 1:n_basefuncs_A
            N_AT = shape_value(cellvalues_A, QPi, i)
            @inbounds for j in 1:n_basefuncs_u
                Bu = tovoigt(shape_symmetric_gradient(cellvalues_u, QPi, j), offdiagscale = 2, order = [1 4 6; 7 2 5; 9 8 3])
                Se[BlockIndex((A▄,u▄), (i,j))] += ctan[2]*(-1.0 * (N_AT ⋅ (eₚ * Bu)) * dΩ)
            end
            @inbounds for j in 1:n_basefuncs_ϕ
                Bgrad = shape_gradient(cellvalues_phi, QPi, j) 
                Se[BlockIndex((A▄,ϕ▄), (i,j))] += ctan[1]*((N_AT ⋅ (κ * Bgrad)) * dΩ) + ctan[2]*((N_AT ⋅ (ϵₜ * Bgrad)) * dΩ)
            end
            @inbounds for j in 1:n_basefuncs_A
                N_A = shape_value(cellvalues_A, QPi, j)
                Se[BlockIndex((A▄,A▄), (i,j))] += ctan[2]*((N_AT ⋅ (κ * N_A)) * dΩ) + ctan[3]*((N_AT ⋅ (ϵₜ * N_A)) * dΩ)
            end
        end

        # calculate residual
        @inbounds for i in 1:n_basefuncs_u
            BuT = tovoigt(shape_symmetric_gradient(cellvalues_u, QPi, i), offdiagscale = 2, order = [1 4 6; 7 2 5; 9 8 3])
            re[BlockIndex((u▄), (i))] -= (BuT ⋅ σ) * dΩ 
        end

        @inbounds for i in 1:n_basefuncs_ϕ
            BgradT = shape_gradient(cellvalues_phi, QPi, i)
            re[BlockIndex((ϕ▄), (i))] -= (BgradT ⋅ D) * dΩ
        end

        @inbounds for i in 1:n_basefuncs_A
            BcurlT = shape_curl(cellvalues_A, QPi, i)
            N_AT = shape_value(cellvalues_A, QPi, i)
            # gauge function: div_Ψ = Bdiv * Bdiv' * A
            div_Ψ = shape_divergence(cellvalues_A, QPi, i) ⋅ function_divergence(cellvalues_A, QPi, A) 
            re[BlockIndex((A▄), (i))] -= (BcurlT ⋅ H - N_AT ⋅ Dp - N_AT ⋅ J + γ * div_Ψ) * dΩ
        end
        
    end # of loop QPi
    
end