# assembly of the macroscopic stiffness matrix
function doassemble_macro(S::SparseMatrixCSC, grid::Grid, dh::DofHandler, mp_b::CorticalBone, mp_m::BoneMarrow, mp_s::SurroundingMediumAir, mp::MicroProblem, mtm::MacroTangentModuli, u::Array{Float64,1}, v::Array{Float64,1}, a::Array{Float64,1}, states::Array{Array{MacroMaterialState{Float64,Array{Float64,1}},1},1}, t::Int, cbnumber_a::Array{Float64,1})
    
    @timeit "parloop" result = pmap(i -> elmt_par!(i, grid, mp_b, mp_m, mp_s, mp, mtm, u, v, a, t, dh, cbnumber_a[i]), 1:length(CellIterator(dh))) 
    
    r = zeros(ndofs(dh))
    assembler = start_assemble(S, r)
    @inbounds for i in 1:length(CellIterator(dh))
        # unpack re, se, eldofs 
        se = result[i][1]
        re = result[i][2]
        eldofs = result[i][3]
        # assembly
        assemble!(assembler, eldofs, re, se)
    end
    
    # states
    @inbounds for (el, state) in enumerate(states) 
        state_el = 0
        @inbounds for i in 1:length(CellIterator(dh)) 
            if(el==result[i][5])
                state_el = i
                break
            end
        end
        @inbounds for j in eachindex(state)
            state[j].σ[:] = result[state_el][4][j].σ[:]
            state[j].ε[:] = result[state_el][4][j].ε[:]
            state[j].E[:] = result[state_el][4][j].E[:]
            state[j].D[:] = result[state_el][4][j].D[:]
            state[j].Dp[:] = result[state_el][4][j].Dp[:]
            state[j].B[:] = result[state_el][4][j].B[:]
            state[j].H[:] = result[state_el][4][j].H[:]
            state[j].J[:] = result[state_el][4][j].J[:]
        end
    end 
    
    result = nothing
    
    return S, r
end;

# macroscopic element routine (input for pmap to parallelize calculations)
function elmt_par!(i::Int, grid::Grid, mp_b::CorticalBone, mp_m::BoneMarrow, mp_s::SurroundingMediumAir, mp::MicroProblem, mtm::MacroTangentModuli, u::Array{Float64,1}, v::Array{Float64,1}, a::Array{Float64,1}, t::Int, dh::DofHandler, cbnumber::Float64)
    # i -> elmtno
    
    celltype = getcelltype(grid)
    cellvalues_u, cellvalues_phi, cellvalues_A = create_values(celltype)
    
    nu = getnbasefunctions(cellvalues_u)
    nϕ = getnbasefunctions(cellvalues_phi)
    nA = getnbasefunctions(cellvalues_A)

    re = PseudoBlockArray(zeros(nu+nϕ+nA), [nu, nϕ, nA]) # local force vector
    se = PseudoBlockArray(zeros(nu+nϕ+nA,nu+nϕ+nA), [nu, nϕ, nA], [nu, nϕ, nA]) # local stiffness matrix
    fill!(se, 0)
    fill!(re, 0)
    
    eldofs = zeros(Int, ndofs_per_cell(dh))
    celldofs!(eldofs, dh, i)

    ue = zeros(ndofs_per_cell(dh))
    ve = zeros(ndofs_per_cell(dh))
    ae = zeros(ndofs_per_cell(dh))
    ue[:] = u[eldofs][:]
    ve[:] = v[eldofs][:]
    ae[:] = a[eldofs][:]
    
    n_qp = getnquadpoints(cellvalues_u)
    state = [MacroMaterialState() for _ in 1:n_qp]

    if(i in getcellset(grid, "air"))
        elmt_air_macro!(se, re, i, dh, cellvalues_phi, cellvalues_A, grid, mp_s, mp.sp, ue, ve, ae, state)
    elseif(i in getcellset(grid, "spongybone"))
        elmt_macro!(se, re, i, dh, cellvalues_u, cellvalues_phi, cellvalues_A, grid, mp_b, mp_m, mp, mtm, ue, ve, ae, state, t, cbnumber)
    else
        print("Error: Element without macro material found")
    end
        
    return se, re, eldofs, state, i
end