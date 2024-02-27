#--- struct for macro tangent quantities
struct MacroTangentModuli{T, S <: AbstractArray{T, 2}} 
    C::S # mechanic stiffness tensor
    ϵₜ::S #  dielectric tensor
    μⁱ::S # inverse of magnetic tensor 
    eₚ::S # piezoelectric tensor
    κ::S # conductivity tensor
end

#--- struct for macro state variables
mutable struct MacroMaterialState{T, S <: AbstractArray{T, 1}}
    σ::S # stress 
    ε::S # strain
    E::S # electric field strength 
    D::S # electric displacement field
    Dp::S # time derivative electric displacement field
    B::S # magnetic flux density
    H::S # magnetic field strength
    J::S # electric current density
end

# constructor for struct MacroMaterialState
function MacroMaterialState()
    return MacroMaterialState(zeros(6),zeros(6),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3),zeros(3))
end

#---  simulation parameters /numerical parameters for time integration
struct SimulationParameters{T, S <: AbstractArray{T, 1}} 
    # time integration parameters
    ρ_∞::T 
    α_f::T
    αₘ::T
    γₐ::T

    Δₜ::T # time increment
    NEWTON_TOL::T # tolerance for convergence check of NR-method
    γ::T # penalty parameter for divergence of magnetic vector potential
    ctan::S # weight vector for time integration
end

# constructor for struct SimulationParameters
function SimulationParameters(ρ_∞::Float64, Δₜ::Float64, NEWTON_TOL::Float64, γ::Float64)
    # time integration - method parameters (Eq.32) - JWH-alpha method Kadapa 2017
    # 0 <= ρ_∞ <= 1, numerical damping increases for small choice of this parameter
    α_f = 1/(1+ρ_∞)
    αₘ = (3-ρ_∞)/(2*(1+ρ_∞))
    γₐ = 0.5 + αₘ - α_f
    ctan = [α_f, αₘ/(γₐ*Δₜ), αₘ^2/(α_f*γₐ^2*Δₜ^2)]
    return SimulationParameters(ρ_∞, α_f, αₘ, γₐ, Δₜ, NEWTON_TOL, γ, ctan)
end