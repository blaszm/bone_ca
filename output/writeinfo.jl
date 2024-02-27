# write file info.txt with simulation information for microscale
function writeinfo_micro(mp_b::CorticalBone, mp_m::BoneMarrow, celltype::DataType, n_dim::Int, n_dofs_n::Int, n_qp::Int, n_el::Int, Ω::Float64, ε_macro::Array{Float64,1}, E_macro::Array{Float64,1}, B_macro::Array{Float64,1}, σ_aver::Array{Float64,1}, D_aver::Array{Float64,1}, Dp_aver::Array{Float64,1}, H_aver::Array{Float64,1}, J_aver::Array{Float64,1}, el_mat1::Int, el_mat2::Int, globgpnumber::Int, a::Float64, b::Float64, n::Int)
        work_dir = pwd()
        cd("Results")
    
        gpnumber = @sprintf("%5.5i", globgpnumber)
        info = string("info_micro_", gpnumber, ".txt")
        io = open(info, "w") 
        println(io, "Simulation information \n")
    
        println(io,"Material parameters \n")
    
        println(io,"Cortical bone")
        println(io,"E_b: $(@sprintf("%.4e", mp_b.E))")
        println(io,"nue_b: $(@sprintf("%.4e", mp_b.ν))")
        println(io,"eps_11_b: $(@sprintf("%.4e", mp_b.ϵ₁))")
        println(io,"eps_33_b: $(@sprintf("%.4e", mp_b.ϵ₃))")
        println(io,"mueinv_b: $(@sprintf("%.4e", mp_b.μᶜ))")
        println(io,"e_p14_b: $(@sprintf("%.4e", mp_b.eₚ14))")
 
        println(io,"\nElastic stiffness tensor")    
        @inbounds for i in 1:length(mp_b.Cᵉ[:,1])
            println(io,"$(@sprintf("%.8e", mp_b.Cᵉ[i,1])) $(@sprintf("%.8e", mp_b.Cᵉ[i,2])) $(@sprintf("%.8e", mp_b.Cᵉ[i,3])) $(@sprintf("%.8e", mp_b.Cᵉ[i,4])) $(@sprintf("%.8e", mp_b.Cᵉ[i,5])) $(@sprintf("%.8e", mp_b.Cᵉ[i,6]))") 
        end
    
        println(io,"\nDielectric tensor")  
        @inbounds for i in 1:length(mp_b.ϵₜ[:,1])
            println(io,"$(@sprintf("%.8e", mp_b.ϵₜ[i,1])) $(@sprintf("%.8e", mp_b.ϵₜ[i,2])) $(@sprintf("%.8e", mp_b.ϵₜ[i,3]))")
        end
    
        println(io,"\nInverse magnetic tensor")
        @inbounds for i in 1:length(mp_b.μⁱ[:,1])
            println(io,"$(@sprintf("%.8e", mp_b.μⁱ[i,1])) $(@sprintf("%.8e", mp_b.μⁱ[i,2])) $(@sprintf("%.8e", mp_b.μⁱ[i,3]))")
        end
    
        println(io,"\nPiezoelectric tensor")       
        @inbounds for i in 1:length(mp_b.eₚ[:,1])
            println(io,"$(@sprintf("%.8e", mp_b.eₚ[i,1])) $(@sprintf("%.8e", mp_b.eₚ[i,2])) $(@sprintf("%.8e", mp_b.eₚ[i,3])) $(@sprintf("%.8e", mp_b.eₚ[i,4])) $(@sprintf("%.8e", mp_b.eₚ[i,5])) $(@sprintf("%.8e", mp_b.eₚ[i,6]))") 
        end
        
        println(io,"\nBone marrow")
        println(io,"E_m: $(@sprintf("%.4e", mp_m.E))")
        println(io,"nue_m: $(@sprintf("%.4e", mp_m.ν))")
        println(io,"eps_11_m: $(@sprintf("%.4e", mp_m.ϵ₁))")
        println(io,"eps_33_m: $(@sprintf("%.4e", mp_m.ϵ₃))")
        println(io,"mueinv_m: $(@sprintf("%.4e", mp_m.μᶜ))")
        println(io,"kappa_1_m: $(@sprintf("%.4e", mp_m.κ₁))")
        println(io,"mue_ve_m: $(@sprintf("%.4e", mp_m.μᵥ))")    
    
        println(io,"\nElastic stiffness tensor")       
        @inbounds for i in 1:length(mp_m.Cᵉ[:,1])
            println(io,"$(@sprintf("%.8e", mp_m.Cᵉ[i,1])) $(@sprintf("%.8e", mp_m.Cᵉ[i,2])) $(@sprintf("%.8e", mp_m.Cᵉ[i,3])) $(@sprintf("%.8e", mp_m.Cᵉ[i,4])) $(@sprintf("%.8e", mp_m.Cᵉ[i,5])) $(@sprintf("%.8e", mp_m.Cᵉ[i,6]))") 
        end
    
        println(io,"\nDielectric tensor")  
        @inbounds for i in 1:length(mp_m.ϵₜ[:,1])
            println(io,"$(@sprintf("%.8e", mp_m.ϵₜ[i,1])) $(@sprintf("%.8e", mp_m.ϵₜ[i,2])) $(@sprintf("%.8e", mp_m.ϵₜ[i,3]))")
        end    
    
        println(io,"\nInverse magnetic tensor")
        @inbounds for i in 1:length(mp_m.μⁱ[:,1])
            println(io,"$(@sprintf("%.8e", mp_m.μⁱ[i,1])) $(@sprintf("%.8e", mp_m.μⁱ[i,2])) $(@sprintf("%.8e", mp_m.μⁱ[i,3]))")
        end
    
        println(io,"\nConductivity tensor")
        @inbounds for i in 1:length(mp_m.κ[:,1])
            println(io,"$(@sprintf("%.8e", mp_m.κ[i,1])) $(@sprintf("%.8e", mp_m.κ[i,2])) $(@sprintf("%.8e", mp_m.κ[i,3]))")
        end
    
        # mesh parameters
        println(io,"\nMesh parameters \n")
        println(io,"Mesh resolution: $n")
        println(io,"Element type: $celltype")
        n_npe = Ferrite.nnodes(celltype) # number of nodes per element
        n_fpe = Ferrite.nfaces(celltype) # number of faces per element
        println(io,"Dim $n_dim, nodes per element $n_npe, faces per element $n_fpe")
        println(io,"DoF per node: $n_dofs_n")
        println(io,"Quadrature points per element: $n_qp")
        println(io,"Total number of elements: $n_el")
        println(io,"Total volume: $(@sprintf("%.4e", Ω))")
    
        println(io,"Number elements of Material 1: $el_mat1")
        println(io,"Number elements of Material 2: $el_mat2")
    
        println(io,"\nRVE type:")
        cbpc = 100*(6*a*b^2+b^3)/((2*a+b)^3)
        println(io,"Vol. perc. of cortical bone: $(@sprintf("%.1f", cbpc))")

        # macro quantities
        println(io,"\nMacro strain:")
        println(io,"$(@sprintf("%.8e", ε_macro[1])) $(@sprintf("%.8e", ε_macro[2])) $(@sprintf("%.8e", ε_macro[3])) $(@sprintf("%.8e", ε_macro[4])) $(@sprintf("%.8e", ε_macro[5])) $(@sprintf("%.8e", ε_macro[6]))")
    
        println(io,"\nMacro electric field:")
        println(io,"$(@sprintf("%.8e", E_macro[1])) $(@sprintf("%.8e", E_macro[2])) $(@sprintf("%.8e", E_macro[3]))")

        println(io,"\nMacro magnetic flux density:")
        println(io,"$(@sprintf("%.8e", B_macro[1])) $(@sprintf("%.8e", B_macro[2])) $(@sprintf("%.8e", B_macro[3]))")
    
        # volume averages
        println(io,"\nVolume average stress:")
        println(io,"$(@sprintf("%.8e", σ_aver[1])) $(@sprintf("%.8e", σ_aver[2])) $(@sprintf("%.8e", σ_aver[3])) $(@sprintf("%.8e", σ_aver[4])) $(@sprintf("%.8e", σ_aver[5])) $(@sprintf("%.8e", σ_aver[6]))")
    
        println(io,"\nVolume average electric displacement:")
        println(io,"$(@sprintf("%.8e", D_aver[1])) $(@sprintf("%.8e", D_aver[2])) $(@sprintf("%.8e", D_aver[3]))")
    
        println(io,"\nVolume average electric displacement time derivative:")
        println(io,"$(@sprintf("%.8e", Dp_aver[1])) $(@sprintf("%.8e", Dp_aver[2])) $(@sprintf("%.8e", Dp_aver[3]))")
    
        println(io,"\nVolume average magnetic field strength:")
        println(io,"$(@sprintf("%.8e", H_aver[1])) $(@sprintf("%.8e", H_aver[2])) $(@sprintf("%.8e", H_aver[3]))")
    
        println(io,"\nVolume average electric current density:")
        println(io,"$(@sprintf("%.8e", J_aver[1])) $(@sprintf("%.8e", J_aver[2])) $(@sprintf("%.8e", J_aver[3]))")

        close(io)
        cd(work_dir)
    return
end

# write file info.txt with simulation information
function writeinfo_macro(mp_b::CorticalBone, mp_m::BoneMarrow, mp_s::SurroundingMediumAir, n_dim::Int, n_npe::Int, n_fpe::Int, n_dofs_n::Int, n_qp::Int, n_el::Int, sp::SimulationParameters, modeltype::String, micromodel::String, celltype_macro::DataType, celltype_micro::DataType)
    
        work_dir = pwd()
        cd("Results")

        io = open("modelinfo.txt", "w") 
        println(io,"Simulation information\n")
    
        println(io,"Material parameters\n")
    
        println(io,"Cortical bone")
        println(io,"E_b: $(@sprintf("%.4e", mp_b.E))")
        println(io,"nue_b: $(@sprintf("%.4e", mp_b.ν))")
        println(io,"eps_11_b: $(@sprintf("%.4e", mp_b.ϵ₁))")
        println(io,"eps_33_b: $(@sprintf("%.4e", mp_b.ϵ₃))")
        println(io,"mueinv_b: $(@sprintf("%.4e", mp_b.μᶜ))")
        println(io,"e_p14_b: $(@sprintf("%.4e", mp_b.eₚ14))")
    
        println(io,"\nElastic stiffness tensor")    
        @inbounds for i in 1:length(mp_b.Cᵉ[:,1])
            println(io,"$(@sprintf("%.8e", mp_b.Cᵉ[i,1])) $(@sprintf("%.8e", mp_b.Cᵉ[i,2])) $(@sprintf("%.8e", mp_b.Cᵉ[i,3])) $(@sprintf("%.8e", mp_b.Cᵉ[i,4])) $(@sprintf("%.8e", mp_b.Cᵉ[i,5])) $(@sprintf("%.8e", mp_b.Cᵉ[i,6]))") 
        end
    
        println(io,"\nDielectric tensor")  
        @inbounds for i in 1:length(mp_b.ϵₜ[:,1])
            println(io,"$(@sprintf("%.8e", mp_b.ϵₜ[i,1])) $(@sprintf("%.8e", mp_b.ϵₜ[i,2])) $(@sprintf("%.8e", mp_b.ϵₜ[i,3]))")
        end

        println(io,"\nInverse magnetic tensor")
        @inbounds for i in 1:length(mp_b.μⁱ[:,1])
            println(io,"$(@sprintf("%.8e", mp_b.μⁱ[i,1])) $(@sprintf("%.8e", mp_b.μⁱ[i,2])) $(@sprintf("%.8e", mp_b.μⁱ[i,3]))")
        end
    
        println(io,"\nPiezoelectric tensor")       
        @inbounds for i in 1:length(mp_b.eₚ[:,1])
            println(io,"$(@sprintf("%.8e", mp_b.eₚ[i,1])) $(@sprintf("%.8e", mp_b.eₚ[i,2])) $(@sprintf("%.8e", mp_b.eₚ[i,3])) $(@sprintf("%.8e", mp_b.eₚ[i,4])) $(@sprintf("%.8e", mp_b.eₚ[i,5])) $(@sprintf("%.8e", mp_b.eₚ[i,6]))") 
        end
        
        println(io,"Bone marrow")
        println(io,"E_m: $(@sprintf("%.4e", mp_m.E))")
        println(io,"nue_m: $(@sprintf("%.4e", mp_m.ν))")
        println(io,"eps_11_m: $(@sprintf("%.4e", mp_m.ϵ₁))")
        println(io,"eps_33_m: $(@sprintf("%.4e", mp_m.ϵ₃))")
        println(io,"mueinv_m: $(@sprintf("%.4e", mp_m.μᶜ))")
        println(io,"kappa_1_m: $(@sprintf("%.4e", mp_m.κ₁))")
        println(io,"mue_ve_m: $(@sprintf("%.4e", mp_m.μᵥ))")    
    
        println(io,"\nElastic stiffness tensor")       
        @inbounds for i in 1:length(mp_m.Cᵉ[:,1])
            println(io,"$(@sprintf("%.8e", mp_m.Cᵉ[i,1])) $(@sprintf("%.8e", mp_m.Cᵉ[i,2])) $(@sprintf("%.8e", mp_m.Cᵉ[i,3])) $(@sprintf("%.8e", mp_m.Cᵉ[i,4])) $(@sprintf("%.8e", mp_m.Cᵉ[i,5])) $(@sprintf("%.8e", mp_m.Cᵉ[i,6]))") 
        end
    
        println(io,"\nDielectric tensor")  
        @inbounds for i in 1:length(mp_m.ϵₜ[:,1])
            println(io,"$(@sprintf("%.8e", mp_m.ϵₜ[i,1])) $(@sprintf("%.8e", mp_m.ϵₜ[i,2])) $(@sprintf("%.8e", mp_m.ϵₜ[i,3]))")
        end    
    
        println(io,"\nInverse magnetic tensor")
        @inbounds for i in 1:length(mp_m.μⁱ[:,1])
            println(io,"$(@sprintf("%.8e", mp_m.μⁱ[i,1])) $(@sprintf("%.8e", mp_m.μⁱ[i,2])) $(@sprintf("%.8e", mp_m.μⁱ[i,3]))")
        end
    
        println(io,"\nConductivity tensor")
        @inbounds for i in 1:length(mp_m.κ[:,1])
            println(io,"$(@sprintf("%.8e", mp_m.κ[i,1])) $(@sprintf("%.8e", mp_m.κ[i,2])) $(@sprintf("%.8e", mp_m.κ[i,3]))")
        end

        println(io,"Surrounding medium air")
        println(io,"eps_11_o: $(@sprintf("%.4e", mp_s.ϵ₁))")
        println(io,"eps_33_o: $(@sprintf("%.4e", mp_s.ϵ₃))")
        println(io,"mueinv_o: $(@sprintf("%.4e", mp_s.μᶜ))")

        println(io,"\nDielectric tensor")  
        @inbounds for i in 1:length(mp_m.ϵₜ[:,1])
            println(io,"$(@sprintf("%.8e", mp_m.ϵₜ[i,1])) $(@sprintf("%.8e", mp_m.ϵₜ[i,2])) $(@sprintf("%.8e", mp_m.ϵₜ[i,3]))")
        end    
    
        println(io,"\nInverse magnetic tensor")
        @inbounds for i in 1:length(mp_m.μⁱ[:,1])
            println(io,"$(@sprintf("%.8e", mp_m.μⁱ[i,1])) $(@sprintf("%.8e", mp_m.μⁱ[i,2])) $(@sprintf("%.8e", mp_m.μⁱ[i,3]))")
        end

        # macro mesh parameters
        println(io,"\nMacro mesh parameters\n")
        println(io,"Element type:")
        println(io,"Dim $n_dim, nodes per element $n_npe, faces per element $n_fpe")
        println(io,"DoF per node: $n_dofs_n")
        println(io,"Quadrature points per element: $n_qp")
        println(io,"Total number of elements: $n_el")
    
        # numerical parameters
        println(io,"\nNumerical parameters\n")
        
        println(io,"Rho_inf: $(@sprintf("%.4e", sp.ρ_∞))")
        println(io,"Delta t: $(@sprintf("%.4e", sp.Δₜ))")
        println(io,"Newton tolerance: $(@sprintf("%.4e", sp.NEWTON_TOL))")
        println(io,"Gauge penalization parameter gamma: $(@sprintf("%.4e", sp.γ))")
    
        # model
        println(io,"\nModel\n")
        println(io,"Model type: $modeltype")
        println(io,"Micromodel: $micromodel")
        println(io,"Element type macro: $celltype_macro")
        println(io,"Element type micro: $celltype_micro")
    
        close(io)
        cd(work_dir)
    return
end