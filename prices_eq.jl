using Pkg
using LinearAlgebra
using Plots
using Distributions
using Parameters, Random, Combinatorics

@with_kw struct parameters
    I::Int64 = 1000      # Number of households
    J::Int64 = 100       # Number of firms
    L::Int64 = 20       # Number of locations
    σ1::Float64 = 1.0   # Utility parameter
    σ2::Float64 = 0.5   # Utility parameter
    ψ::Float64 = 0.25   # Work distaste parameter
    τ1::Float64 = 0     # Fixed cost of shopping
    τ2::Float64 = 0.5   # Variable cost of shopping
    γ1::Float64 = 0     # Fixed cost of commuting
    γ2::Float64 = 0.5   # Variable cost of commuting
    κ::Float64 = 5      # Non-sunk firm fixed cost
    α::Float64 = 0.7    # Production function DRTS parameter
    φ::Float64 = 0.4    # Wage adjustment increment
    u0::Float64 = 0     # Outside option utility
    λ::Float64 = 0.5    # Marginal cost of land
    ρ::Float64 = 3      # Shape parameter for Pareto prod. multipliers
    θ0::Float64 = 0.25  # Minimum prod. multiplier
end

para = parameters()

function price_eq(para, Kl, hhs, hl, livein, firms, active, offers_guess, prices_guess, Θ, it_tols)
    @unpack J, L, ψ, κ, α = para
    # Initialize while loop
    iter = 0
    maxdiff = 1e+6
    prices = prices_guess
    prices_store = zeros(J,Int(it_tols[3])+2)
    prices_store[:,2] = prices
    offers = offers_guess
    Z = zeros(size(active)[1])
    matches = zeros(size(active)[1],size(livein)[1])
    Ω_star = []
    while (maxdiff > it_tols[4]) & (iter < it_tols[3])
        iter += 1
        offers_guess = zeros(J,I,2)
        offers_guess[:,:,1] .= 10
        offers_guess[:,:,2] .= 1
        # Find wage equilibrium given current guess
        offers, Z, matches, Ω_star, Ωhj = wage_eq(para, Kl, hhs, hl, livein, firms, active, offers_guess, prices, Θ, it_tols)

        # Calculate marginal costs
        prices_new = zeros(size(active)[1])
        for j in active
           if size(Ω_star[j])[1] > 0
               prices_new[j] = ((Z[j]/Kl[firms[j]])^((1-α)/α))*(1/(α*Kl[firms[j]]*sum(Θ[j,Ω_star[j]])))*transpose(offers[j,Ω_star[j],1])*offers[j,Ω_star[j],2] + 1/σ1
           else
               prices_new[j] = prices[j]
           end #size check
        end # j-loop

        prices_store[:,iter+2] = prices_new
        # calculate maximum difference between prices
        maxdiff = findmax(abs.(prices - prices_new))[1]

        # Check for convergence
        if maxdiff < it_tols[4]
            println("Price equilibrium achieved after $iter iterations")
        else
            prices = 0.4*(prices_new + prices) +
            0.2*(prices_store[:,iter + 1] + prices_store[:,iter])
            println("Price Iteration $iter, $maxdiff")
        end # convergence check
    end # price convergence loop
    return offers, Z, matches, Ω_star, prices
end # price_eq

prices_guess = ones(J)
offers_ret, Z_ret, matches_ret, Ostar, prices_ret = price_eq(para, Kl, hhs, hl, livein, firms, active, offers_guess, prices_guess, Θ, it_tols)
