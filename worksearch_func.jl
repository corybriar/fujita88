using Pkg
using LinearAlgebra
using Plots
using Distributions
using Parameters, Random

@with_kw struct parameters
    I::Int64 = 5000     # Number of agents
    J::Int64 = 100        # Number of firms
    L::Int64 = 20       # Number of locations
    σ1::Float64 = 0.5   # Utility parameter
    σ2::Float64 = 0.5  # Utility parameter
    ψ::Float64 = 0.25   # Work distaste parameter
    τ1::Float64 = 0.5   # Fixed cost of shopping
    τ2::Float64 = 0.5   # Variable cost of shopping
    γ1::Float64 = 0.5   # Fixed cost of commuting
    γ2::Float64 = 0.5   # Variable cost of commuting
    κ::Float64 = 5      # Non-sunk firm fixed cost
    α::Float64 = 0.7    # Production function DRTS parameter
    u0::Float64 = 0     # Outside option utility
end

para = parameters()

function τ(para, firms)
    @unpack J, L, τ1, τ2 = para
    locs = zeros(J,L)
    # matrix of locations
    for l in 1:L
        locs[:,l] .= l
    end # l-loop
    # Create matrix of travel costs
    tcosts = τ1 .+ τ2.*abs.(locs .- transpose(firms))
    return tcosts
end #τ


function wage_eq(para, rents, hhs, hl, firms, fl, offers_guess, prices, maxit_wage, wage_tol)
    @unpack I, J, L, σ1, σ2, ψ, τ1, τ2, γ1, γ2, κ, α = para
    # Calculate firm demands
    Z = σ1*σ2.*exp.(-σ1.*(transpose(prices) .+ τ(para,firms)))*transpose(hl)

    # initialize while loop
    iter = 0
    maxdiff = 1e+6
    offers = offers_guess

    while (maxdiff > wage_tol) & (iter < maxit_wage + 1)
        iter += 1

        # Households rank offers
        

    end # convergence loop
    return offers, Z
end # wage_eq


"""
MATRICES:
rents[1,L] -- rents for each location l
hhs[1,I] -- locations of each household (0 = outside)
firms[1,J] -- locations of each firm (0 = not operating)
offers[J,I,2] -- offers by firm j to household i in form (w,n)
prices[1,J] -- prices of goods by each firm j
hl[1,L+1] -- distribution of hh's (L+1 = outside)
fl[1,L+1] -- distribution of firms (L+1 = not operating)

"""