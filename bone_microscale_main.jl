# function to perform calculation considering the microscale only

# showcase example for 2 different RVEs (healthy and ill bone)
# calculates effective stiffness and solves microscale problem 
# for small shear deformation for both cases

const n_timesteps = 1
include("_include.jl")

function microscale_solver()
    reset_timer!()
    
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
    
        return mp_b, mp_m
    end

    # preperation export
    if(isdir("Results") == false)
        mkdir("Results")
    end

    sp = SimulationParameters(0.1, 1e-4, 1e-8, 1.0) # ρ_∞, Δₜ, NEWTON_TOL, γ
    Δₜ = sp.Δₜ

    # micro material
    mp_b, mp_m  = initmaterialparameters(Δₜ)
    
    #------------------------------
    # 1st calculation: healthy bone
    # micro mesh parameters
    a = 0.00032 
    b = 0.00036
    n = 2

    # create micro problem
    mp = MicroProblem(a, b, n, sp, mp_b, mp_m, "rve30pc") 
    
    # preperation HDF files
    hdfname_i = "epsidata_i.h5"
    hdfname_ii = "epsidata_ii.h5"
    n_cells = length(mp.ctx.grid.cells) 
    preparehdf(hdfname_i, hdfname_ii, 8*n_cells, 1, 1)
        
    # calculate effective stiffness from macro tangent moduli
    mtm = calculatemacrotangent(mp, mp_b, mp_m)
    C_12 = mtm.C[1,2]
    C_44 = mtm.C[4,4]
    Em = (C_44*(3*C_12+2*C_44)) / (C_12+C_44)
    println("\nYoungs modulus RVE from numerical perturbation - 30pc bone: $(@sprintf("%.4e", Em))")
    
    # microscale calculation
    ε_macro = zeros(6)
    E_macro = zeros(3)
    B_macro = zeros(3)    
    ε_macro[5] = 0.0002
    mq = MacroQuantities(ε_macro, E_macro, B_macro, 1)
    σ, D, Dp, H, J = solve_RVE(mq, mp, mp_b, mp_m, 1; overwrite=false, output=true)

    println("\nVolume average stress σ:")
    println("$(@sprintf("%.4e", σ[1])) $(@sprintf("%.4e", σ[2])) $(@sprintf("%.4e", σ[3])) $(@sprintf("%.4e", σ[4])) $(@sprintf("%.4e", σ[5])) $(@sprintf("%.4e", σ[6])) \n")

    println("Volume average electric displacement D:")
    println("$(@sprintf("%.4e", D[1])) $(@sprintf("%.4e", D[2])) $(@sprintf("%.4e", D[3])) \n")

    println("Volume average electric displacement time derivative Dp:")
    println("$(@sprintf("%.4e", Dp[1])) $(@sprintf("%.4e", Dp[2])) $(@sprintf("%.4e", Dp[3])) \n")

    println("Volume average magnetic field strength H:")
    println("$(@sprintf("%.4e", H[1])) $(@sprintf("%.4e", H[2])) $(@sprintf("%.4e", H[3])) \n")

    println("Volume average electric current density J:")
    println("$(@sprintf("%.4e", J[1])) $(@sprintf("%.4e", J[2])) $(@sprintf("%.4e", J[3])) \n")


    #------------------------------
    # 2nd calculation: ill bone
    # micro mesh parameters
    a = 0.00043 
    b = 0.00014
    n = 2

    # create micro problem
    mp = MicroProblem(a, b, n, sp, mp_b, mp_m, "rve5pc") 
   
    # preperation HDF files
    hdfname_i = "epsidata_i.h5"
    hdfname_ii = "epsidata_ii.h5"
    n_cells = length(mp.ctx.grid.cells) 
    preparehdf(hdfname_i, hdfname_ii, 8*n_cells, 1, 1)
       
    # calculate effective stiffness from macro tangent moduli
    mtm = calculatemacrotangent(mp, mp_b, mp_m)
    C_12 = mtm.C[1,2]
    C_44 = mtm.C[4,4]
    Em = (C_44*(3*C_12+2*C_44)) / (C_12+C_44)
    println("\nYoungs modulus RVE from numerical perturbation - 5pc bone: $(@sprintf("%.4e", Em))")
   
    # microscale calculation
    ε_macro = zeros(6)
    E_macro = zeros(3)
    B_macro = zeros(3)    
    ε_macro[5] = 0.0002
    mq = MacroQuantities(ε_macro, E_macro, B_macro, 1)
    σ, D, Dp, H, J = solve_RVE(mq, mp, mp_b, mp_m, 1; overwrite=false, output=true)

    println("\nVolume average stress σ:")
    println("$(@sprintf("%.4e", σ[1])) $(@sprintf("%.4e", σ[2])) $(@sprintf("%.4e", σ[3])) $(@sprintf("%.4e", σ[4])) $(@sprintf("%.4e", σ[5])) $(@sprintf("%.4e", σ[6])) \n")

    println("Volume average electric displacement D:")
    println("$(@sprintf("%.4e", D[1])) $(@sprintf("%.4e", D[2])) $(@sprintf("%.4e", D[3])) \n")

    println("Volume average electric displacement time derivative Dp:")
    println("$(@sprintf("%.4e", Dp[1])) $(@sprintf("%.4e", Dp[2])) $(@sprintf("%.4e", Dp[3])) \n")

    println("Volume average magnetic field strength H:")
    println("$(@sprintf("%.4e", H[1])) $(@sprintf("%.4e", H[2])) $(@sprintf("%.4e", H[3])) \n")

    println("Volume average electric current density J:")
    println("$(@sprintf("%.4e", J[1])) $(@sprintf("%.4e", J[2])) $(@sprintf("%.4e", J[3])) \n")
end

microscale_solver()