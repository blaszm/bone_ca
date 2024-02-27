# calculation of elementwise residual and stiffness matrix for cortical bone material model

# ue = element displacement solution vector (at nodes)
# ve = element velocity solution vector (at nodes)
# ae = element acceleration solution vector (at nodes)
# ue stores for each node first all displacements u, then electric scalar potential phi, then magnetic vector potential A

function elmt_cb_r!(re::PseudoBlockArray{Float64,1}, elmt::Ferrite.CellCache, dh::DofHandler, cellvalues_u::CellVectorValues{dim}, cellvalues_phi::CellScalarValues{dim}, cellvalues_A::CellVectorValues{dim}, mp::CorticalBone, sp::SimulationParameters, ue::Array{Float64,1}, ve::Array{Float64,1}, ae::Array{Float64,1}, state::Array{MaterialState{Float64,Array{Float64,1}},1}, mq::MacroQuantities) where {dim}
    n_basefuncs_u = getnbasefunctions(cellvalues_u) # number of shape functions times dimension
    n_basefuncs_ϕ = getnbasefunctions(cellvalues_phi)
    n_basefuncs_A = getnbasefunctions(cellvalues_A)
    
    n_qp = getnquadpoints(cellvalues_u) # number of quadrature points
    
    u▄, ϕ▄ ,A▄ = 1,2,3 # block index for coupled problems
    reinit!(cellvalues_u, elmt)
    reinit!(cellvalues_phi, elmt)
    reinit!(cellvalues_A, elmt)
    
    # unpack simulation parameters
    ctan = zeros(3)
    ε_macro = zeros(6)
    E_macro = zeros(3)
    B_macro = zeros(3)
    
    ctan[:] = sp.ctan[:]
    γ = sp.γ
    
    ε_macro[:] = mq.ε_macro[:]
    E_macro[:] = mq.E_macro[:]
    B_macro[:] = mq.B_macro[:]
    
    # preallocate / reset helper variables
    Cᵉ = mp.Cᵉ 
    ϵₜ = mp.ϵₜ
    μⁱ = mp.μⁱ
    eₚ = mp.eₚ
    
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
        ε += ε_macro
        εp = tovoigt(symmetric(function_gradient(cellvalues_u, QPi, up)), offdiagscale = 2, order = [1 4 6; 7 2 5; 9 8 3])
        
        # electric field
        E = -1.0*function_gradient(cellvalues_phi, QPi, ϕ) - 1.0*function_value(cellvalues_A, QPi, Ap)
        E += E_macro
        
        # magnetic flux density
        B = function_curl(cellvalues_A, QPi, A)
        B += B_macro
    
        # stress 
        σ = Cᵉ * ε - eₚ' * E
        
        # electric displacement field + time derivative
        D = ϵₜ * E + eₚ * ε
        Dp = ϵₜ * (-1.0*function_gradient(cellvalues_phi, QPi, ϕp) - 1.0*function_value(cellvalues_A, QPi, App)) + eₚ * εp
        
        # magnetic field strength
        H = μⁱ * B

        # save material state
        state[QPi].temp_σ = σ
        state[QPi].temp_εⁱ = zeros(6) 
        state[QPi].ε = ε
        state[QPi].E = E
        state[QPi].D = D
        state[QPi].Dp = Dp
        state[QPi].B = B
        state[QPi].H = H
        state[QPi].J = zeros(3)

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
            re[BlockIndex((A▄), (i))] -= (BcurlT ⋅ H - N_AT ⋅ Dp + γ * div_Ψ) * dΩ
        end

    end # of loop QPi
    
end

function elmt_cb_S!(Se::PseudoBlockArray{Float64,2}, elmt::Ferrite.CellCache, cellvalues_u::CellVectorValues{dim}, cellvalues_phi::CellScalarValues{dim}, cellvalues_A::CellVectorValues{dim}, mp::CorticalBone, sp::SimulationParameters) where {dim}   
    n_basefuncs_u = getnbasefunctions(cellvalues_u) # number of shape functions
    n_basefuncs_ϕ = getnbasefunctions(cellvalues_phi)
    n_basefuncs_A = getnbasefunctions(cellvalues_A)
    
    n_qp = getnquadpoints(cellvalues_u) # number of quadrature points
    
    u▄, ϕ▄ ,A▄ = 1,2,3 # block index for coupled problems
    reinit!(cellvalues_u, elmt)
    reinit!(cellvalues_phi, elmt)
    reinit!(cellvalues_A, elmt)
    
    # unpack simulation parameters
    ctan = zeros(3)
    ctan[:] = sp.ctan[:]
    γ = sp.γ
    
    # preallocate / reset helper variables
    Cᵉ = mp.Cᵉ 
    ϵₜ = mp.ϵₜ
    μⁱ = mp.μⁱ
    eₚ = mp.eₚ
    
    @inbounds for QPi in 1:n_qp # loop QPi
        
        # volume part
        dΩ = getdetJdV(cellvalues_u, QPi)
        
        # calculate stiffness / damping / mass matrices

        # Ke_uu += (BuT * Cᵉ * Bu) * dΩ
        # Ke_uϕ += (BuT * eₚ' * Bgrad) * dΩ
        # Ce_uA += (BuT * eₚ' * N_A) * dΩ
        @inbounds for i in 1:n_basefuncs_u
            BuT = tovoigt(shape_symmetric_gradient(cellvalues_u, QPi, i), offdiagscale = 2, order = [1 4 6; 7 2 5; 9 8 3])'
            @inbounds for j in 1:n_basefuncs_u
                Bu = tovoigt(shape_symmetric_gradient(cellvalues_u, QPi, j), offdiagscale = 2, order = [1 4 6; 7 2 5; 9 8 3])
                Se[BlockIndex((u▄,u▄), (i,j))] += ctan[1]*((BuT * Cᵉ * Bu) * dΩ)
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
        # Ce_Aϕ += (N_AT ⋅ (ϵₜ * Bgrad)) * dΩ
        # Me_AA += (N_AT ⋅ (ϵₜ * N_A)) * dΩ
        @inbounds for i in 1:n_basefuncs_A
            N_AT = shape_value(cellvalues_A, QPi, i)
            @inbounds for j in 1:n_basefuncs_u
                Bu = tovoigt(shape_symmetric_gradient(cellvalues_u, QPi, j), offdiagscale = 2, order = [1 4 6; 7 2 5; 9 8 3])
                Se[BlockIndex((A▄,u▄), (i,j))] += ctan[2]*(-1.0 * (N_AT ⋅ (eₚ * Bu)) * dΩ)
            end
            @inbounds for j in 1:n_basefuncs_ϕ
                Bgrad = shape_gradient(cellvalues_phi, QPi, j) 
                Se[BlockIndex((A▄,ϕ▄), (i,j))] += ctan[2]*((N_AT ⋅ (ϵₜ * Bgrad)) * dΩ)
            end
            @inbounds for j in 1:n_basefuncs_A
                N_A = shape_value(cellvalues_A, QPi, j)
                Se[BlockIndex((A▄,A▄), (i,j))] += ctan[3]*((N_AT ⋅ (ϵₜ * N_A)) * dΩ)
            end
        end
        
    end # of loop QPi

end