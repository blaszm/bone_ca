# assembly of the microscale system
# note: Material_1 = cortical bone, Material_2 = bone marrow
function doassemble_r(cellvalues_u::CellVectorValues{dim}, cellvalues_phi::CellScalarValues{dim}, cellvalues_A::CellVectorValues{dim}, grid::Grid, dh::DofHandler, mp_b::CorticalBone, mp_m::BoneMarrow, sp::SimulationParameters, u::Array{Float64,1}, v::Array{Float64,1}, a::Array{Float64,1}, states::Array{Array{MaterialState{Float64,Array{Float64,1}},1},1}, mq::MacroQuantities) where {dim}

    r = zeros(ndofs(dh))
    nu = getnbasefunctions(cellvalues_u)
    nϕ = getnbasefunctions(cellvalues_phi)
    nA = getnbasefunctions(cellvalues_A)
    re = PseudoBlockArray(zeros(nu+nϕ+nA), [nu, nϕ, nA]) # local residual vector

    @inbounds for (elmt, state) in zip(CellIterator(dh), states)

        fill!(re, 0)
        eldofs = celldofs(elmt)
        ue = u[eldofs]
        ve = v[eldofs]
        ae = a[eldofs]
        
        if cellid(elmt) in getcellset(grid, "Material_1")
            elmt_cb_r!(re, elmt, dh, cellvalues_u, cellvalues_phi, cellvalues_A, mp_b, sp, ue, ve, ae, state, mq)
        elseif cellid(elmt) in getcellset(grid, "Material_2")
            elmt_bm_r!(re, elmt, dh, cellvalues_u, cellvalues_phi, cellvalues_A, mp_m, sp, ue, ve, ae, state, mq)
        else
            print("Error: Element without Material found")
        end
        
        assemble!(r, eldofs, re)
    end
    return r
end

function doassemble_S!(cellvalues_u::CellVectorValues{dim}, cellvalues_phi::CellScalarValues{dim}, cellvalues_A::CellVectorValues{dim},
                    S::SparseMatrixCSC, grid::Grid, dh::DofHandler, mp_b::CorticalBone, mp_m::BoneMarrow, sp::SimulationParameters) where {dim}
    
    assembler = start_assemble(S)
    nu = getnbasefunctions(cellvalues_u)
    nϕ = getnbasefunctions(cellvalues_phi)
    nA = getnbasefunctions(cellvalues_A)
    Se = PseudoBlockArray(zeros(nu+nϕ+nA,nu+nϕ+nA), [nu, nϕ, nA], [nu, nϕ, nA]) # local stiffness matrix
    
    @inbounds for elmt in CellIterator(dh)
        fill!(Se, 0)
        eldofs = celldofs(elmt)
        if cellid(elmt) in getcellset(grid, "Material_1")
            elmt_cb_S!(Se, elmt, cellvalues_u, cellvalues_phi, cellvalues_A, mp_b, sp)
        elseif cellid(elmt) in getcellset(grid, "Material_2")
            elmt_bm_S!(Se, elmt, cellvalues_u, cellvalues_phi, cellvalues_A, mp_m, sp)
        else
            print("Error: Element without Material found")
        end
        
        assemble!(assembler, eldofs, Se)
    end

    return
end