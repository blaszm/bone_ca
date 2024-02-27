# --- cortical bone
struct CorticalBone{T, S <: AbstractArray{T, 2}} # T = real (number type)
    E::T # Youngs modulus
    ν::T # Poisson ratio
    Cᵉ::S # elastic stiffness tensor
    ϵ₁::T # dielectric constant 11
    ϵ₃::T # dielectric constant 33
    ϵₜ::S #  dielectric tensor
    μᶜ::T # inverse of magnetic constant
    μⁱ::S # inverse of magnetic tensor 
    eₚ14::T # dielectric constant 14
    eₚ::S # piezoelectric tensor
end

# constructor for struct CorticalBone
function CorticalBone(E::Float64, ν::Float64, ϵ₁::Float64, ϵ₃::Float64, μᶜ::Float64, eₚ14::Float64) 
    # elastic tensor
    Cᵉ = zeros(6,6)
    
    Cᵉ[1,1] = 1 - ν
    Cᵉ[1,2] = ν
    Cᵉ[1,3] = ν
    
    Cᵉ[2,1] = ν
    Cᵉ[2,2] = 1 - ν
    Cᵉ[2,3] = ν
    
    Cᵉ[3,1] = ν
    Cᵉ[3,2] = ν
    Cᵉ[3,3] = 1 - ν
    
    Cᵉ[4,4] = (1 - 2ν)/2
    
    Cᵉ[5,5] = (1 - 2ν)/2
    
    Cᵉ[6,6] = (1 - 2ν)/2
    
    Cᵉ = E/((1 + ν)*(1 - 2ν)) * Cᵉ
    
    # dielectric tensor
    ϵₜ = zeros(3,3)
    
    ϵₜ[1,1] = ϵ₁
    ϵₜ[2,2] = ϵ₁
    ϵₜ[3,3] = ϵ₃
    
    # magnetic tensor
    μⁱ = zeros(3,3)
    
    μⁱ[1,1] = μᶜ
    μⁱ[2,2] = μᶜ
    μⁱ[3,3] = μᶜ
    
    # piezoelectric tensor
    eₚ = zeros(3,6)
    
    eₚ[1,5] = eₚ14
    eₚ[2,6] = -eₚ14 
    
    return CorticalBone(E, ν, Cᵉ, ϵ₁, ϵ₃, ϵₜ, μᶜ, μⁱ, eₚ14, eₚ)
end


# --- bone marrow
struct BoneMarrow{T, S <: AbstractArray{T, 2}} 
    E::T # Youngs modulus
    ν::T # Poisson ratio
    Cᵉ::S # Elastic stiffness tensor
    ϵ₁::T # dielectric constant 11
    ϵ₃::T # dielectric constant 33
    ϵₜ::S #  dielectric tensor
    μᶜ::T # inverse of magnetic constant
    μⁱ::S # inverse of magnetic tensor
    κ₁::T # conductivity constant
    κ::S # conductivity tensor
    μᵥ::T # viscoelastic damping parameter (=delta_t*r_1)
end

# constructor for struct BoneMarrow
function BoneMarrow(E::Float64, ν::Float64, ϵ₁::Float64, ϵ₃::Float64, μᶜ::Float64, κ₁::Float64, μᵥ::Float64) 
    # elastic tensor
    Cᵉ = zeros(6,6)
    
    Cᵉ[1,1] = 1 - ν
    Cᵉ[1,2] = ν
    Cᵉ[1,3] = ν
    
    Cᵉ[2,1] = ν
    Cᵉ[2,2] = 1 - ν
    Cᵉ[2,3] = ν
    
    Cᵉ[3,1] = ν
    Cᵉ[3,2] = ν
    Cᵉ[3,3] = 1 - ν
    
    Cᵉ[4,4] = (1 - 2ν)/2
    
    Cᵉ[5,5] = (1 - 2ν)/2
    
    Cᵉ[6,6] = (1 - 2ν)/2
    
    Cᵉ = E/((1 + ν)*(1 - 2ν)) * Cᵉ
    
    # dielectric tensor
    ϵₜ = zeros(3,3)
    
    ϵₜ[1,1] = ϵ₁
    ϵₜ[2,2] = ϵ₁
    ϵₜ[3,3] = ϵ₃
    
    # magnetic tensor
    μⁱ = zeros(3,3)
    
    μⁱ[1,1] = μᶜ
    μⁱ[2,2] = μᶜ
    μⁱ[3,3] = μᶜ
    
    # conductivity tensor
    κ = zeros(3,3)
    
    κ[1,1] = κ₁
    κ[2,2] = κ₁
    κ[3,3] = κ₁
    
    return BoneMarrow(E, ν, Cᵉ, ϵ₁, ϵ₃, ϵₜ, μᶜ, μⁱ, κ₁, κ, μᵥ)
end

# --- surrounding medium air
struct SurroundingMediumAir{T, S <: AbstractArray{T, 2}} 
    ϵ₁::T # dielectric constant 11
    ϵ₃::T # dielectric constant 33
    ϵₜ::S #  dielectric tensor
    μᶜ::T # inverse of magnetic constant
    μⁱ::S # inverse of magnetic tensor
end

# constructor for struct SurroundingMediumAir
function SurroundingMediumAir(ϵ₁::Float64, ϵ₃::Float64, μᶜ::Float64) 
    
    # dielectric tensor
    ϵₜ = zeros(3,3)
    
    ϵₜ[1,1] = ϵ₁
    ϵₜ[2,2] = ϵ₁
    ϵₜ[3,3] = ϵ₃
    
    # magnetic tensor
    μⁱ = zeros(3,3)
    
    μⁱ[1,1] = μᶜ
    μⁱ[2,2] = μᶜ
    μⁱ[3,3] = μᶜ
    
    return SurroundingMediumAir(ϵ₁, ϵ₃, ϵₜ, μᶜ, μⁱ)
end

# --- state variables
mutable struct MaterialState{T, S <: AbstractArray{T, 1}}
    # store "converged" values of stress and inelastic strain
    σ::S # stress
    εⁱ::S # inelastic strain
    
    # store temporary values used during equilibrium iterations (only necessary for σ, εⁱ)
    temp_σ::S
    temp_εⁱ::S
    
    # other state variables
    ε::S # strain
    
    E::S # electric field strength 
    D::S # electric displacement field
    Dp::S # time derivative electric displacement field
    
    B::S # magnetic flux density
    H::S # magnetic field strength
    
    J::S # electric current density
end

# constructor for struct MaterialState
function MaterialState() 
    return MaterialState(zeros(6),zeros(6),zeros(6),zeros(6),zeros(6),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3))
end

# update function
function update_state!(state::MaterialState)
    state.σ = state.temp_σ
    state.εⁱ = state.temp_εⁱ
end

# struct macroscopic quantities 
struct MacroQuantities{T, S <: AbstractArray{T, 1}}
    ε_macro::S # macro strain
    E_macro::S # macro electric field
    B_macro::S # macro magnetic flux density
    t::Int64 # macro time
end