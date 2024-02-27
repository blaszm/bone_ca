# functions to collect quantities from material states for visualization

function collect_scalar_quantity(dh::DofHandler, states, symbol, i::Int) 
    
    # i = tensor component, order: x,y,z,xy,yz,xz - same as PV default
    field_out = [Vec{1,Float64}[] for _ in 1:getncells(dh.grid)] # Array(Array(Tensor))
    field_i = [0.0] # one component array

    @inbounds for (elmt, cell_states) in enumerate(states)
        field_cell = field_out[elmt]
        
        for state in cell_states
            field_i[1] = getfield(state,symbol)[i]
            field_t = reinterpret(Vec{1, Float64}, vec(field_i)) # convert array to tensor
            push!(field_cell, field_t[1]) # add tensor to Array(Array()) 
        end
    end

    return field_out
end

function collect_vector_quantity(dh::DofHandler, states, symbol) 
    field_out = [Vec{3,Float64}[] for _ in 1:getncells(dh.grid)] 
    field_i = [0.0, 0.0, 0.0]

    @inbounds for (elmt, cell_states) in enumerate(states)
        field_cell = field_out[elmt]
        
        @inbounds for state in cell_states
            field_i[:] = getfield(state,symbol)
            field_t = reinterpret(Vec{3, Float64}, vec(field_i[:]))
            push!(field_cell, field_t[1]) 
        end
    end

    return field_out
end

function collect_tensor_quantity(dh::DofHandler, states, symbol) 
    field_out = [Tensor{2, 3, Float64, 9}[] for _ in 1:getncells(dh.grid)] 

    @inbounds for (elmt, cell_states) in enumerate(states)
        field_cell = field_out[elmt]
        
        for state in cell_states
            
            # convert voigt to tensor notation
            # i = tensor component, order: x,y,z,xy,yz,xz - same as PV default
            field_voigt = getfield(state,symbol)
            field_array = [field_voigt[1] field_voigt[4] field_voigt[6]; field_voigt[4] field_voigt[2] field_voigt[5]; field_voigt[6] field_voigt[5] field_voigt[3]]

            # type conversion to tensor
            field = reinterpret(Tensor{2, 3, Float64, 9}, vec(field_array))

            push!(field_cell, field[1]) 
        end
    end

    return field_out
end