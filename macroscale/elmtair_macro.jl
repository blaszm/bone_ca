# calculation of elementwise residual and stiffness matrix for (macroscale surrounding medium) air material model

# ue = element displacement solution vector (at nodes)
# ve = element velocity solution vector (at nodes)
# ae = element acceleration solution vector (at nodes)
# ue stores for each node first all displacements u, then electric scalar potential phi, then magnetic vector potential A

function elmt_air_macro!(Se::PseudoBlockArray{Float64,2}, re::PseudoBlockArray{Float64,1}, i::Int, dh::DofHandler, cellvalues_phi::CellScalarValues{dim}, cellvalues_A::CellVectorValues{dim}, grid::Grid, mp::SurroundingMediumAir, sp::SimulationParameters, ue::Array{Float64,1}, ve::Array{Float64,1}, ae::Array{Float64,1}, state::Array{MacroMaterialState{Float64,Array{Float64,1}},1}) where {dim} 
    n_basefuncs_ϕ = getnbasefunctions(cellvalues_phi) # number of shape functions times dimension
    n_basefuncs_A = getnbasefunctions(cellvalues_A)
    
    n_qp = getnquadpoints(cellvalues_A) # number of quadrature points
    
    coordinates = [zero(Vec{3}) for _ in 1:length(grid.cells[1].nodes)]
    
    nodeids = grid.cells[i].nodes
    @inbounds for j in eachindex(coordinates)
        coordinates[j] = grid.nodes[nodeids[j]].x
    end 
    
    ϕ▄ ,A▄ = 2,3 # block index for coupled problems
    reinit!(cellvalues_phi, coordinates)
    reinit!(cellvalues_A, coordinates)
    
    # unpack simulation parameters
    ctan = zeros(3)
    ctan[:] = sp.ctan[:]
    γ = sp.γ
       
    # preallocate / reset helper variables
    ϵₜ = mp.ϵₜ
    μⁱ = mp.μⁱ
    
    # main variables (nodes)
    ϕ = zeros(n_basefuncs_ϕ)
    A = zeros(n_basefuncs_A)
    ϕp = zeros(n_basefuncs_ϕ) 
    Ap = zeros(n_basefuncs_A)
    App = zeros(n_basefuncs_A)
    
    ϕ[:] = ue[dof_range(dh, :phi)] 
    A[:] = ue[dof_range(dh, :A)]
    ϕp[:] = ve[dof_range(dh, :phi)] 
    Ap[:] = ve[dof_range(dh, :A)] 
    App[:] = ae[dof_range(dh, :A)] 
    
    @inbounds for QPi in 1:n_qp # loop QPi

        # volume part
        dΩ = getdetJdV(cellvalues_A, QPi)
        
        # electric field
        E = -1.0*function_gradient(cellvalues_phi, QPi, ϕ) - 1.0*function_value(cellvalues_A, QPi, Ap)
        
        # magnetic flux density
        B = function_curl(cellvalues_A, QPi, A)
        
        # electric displacement field + time derivative
        D = ϵₜ * E
        Dp = ϵₜ * (-1.0*function_gradient(cellvalues_phi, QPi, ϕp) - 1.0*function_value(cellvalues_A, QPi, App))
        
        # magnetic field strength
        H = μⁱ * B
        
        # save material state
        state[QPi].σ = zeros(6)
        state[QPi].ε = zeros(6)
        state[QPi].E = E
        state[QPi].D = D
        state[QPi].Dp = Dp
        state[QPi].B = B
        state[QPi].H = H
        state[QPi].J = zeros(3)

        # calculate stiffness / damping / mass matrices

        # Ke_ϕϕ += -1.0 * (BgradT ⋅ (ϵₜ * Bgrad)) * dΩ  
        # Ce_ϕA += -1.0 * (BgradT ⋅ (ϵₜ * N_A)) * dΩ
        @inbounds for i in 1:n_basefuncs_ϕ
            BgradT = shape_gradient(cellvalues_phi, QPi, i)
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

        # Ce_Aϕ += (N_AT ⋅ (ϵₜ * Bgrad)) * dΩ
        # Me_AA += (N_AT ⋅ (ϵₜ * N_A)) * dΩ
        @inbounds for i in 1:n_basefuncs_A
            N_AT = shape_value(cellvalues_A, QPi, i)
            @inbounds for j in 1:n_basefuncs_ϕ
                Bgrad = shape_gradient(cellvalues_phi, QPi, j) 
                Se[BlockIndex((A▄,ϕ▄), (i,j))] += ctan[2]*((N_AT ⋅ (ϵₜ * Bgrad)) * dΩ)
            end
            @inbounds for j in 1:n_basefuncs_A
                N_A = shape_value(cellvalues_A, QPi, j)
                Se[BlockIndex((A▄,A▄), (i,j))] += ctan[3]*((N_AT ⋅ (ϵₜ * N_A)) * dΩ)
            end
        end
 
        # calculate residual
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