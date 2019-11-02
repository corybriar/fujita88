using Pkg
using LinearAlgebra
using Plots
using Distributions
using Parameters, Random

@with_kw struct parameters
    I::Int64 = 5000     # Number of agents
    J::Int64 = 100      # Number of firms
    L::Int64 = 20       # Number of locations
    σ1::Float64 = 0.5   # Utility parameter
    σ2::Float64 = 0.5   # Utility parameter
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

function price_eq(para, Kl, rents, hhs, hl, firms, fl, offers_guess, prices_guess, order, it_tols)
    @unpack J, L, ψ, κ, α = para
    # Initialize while loop
    iter = 0
    maxdiff = 1e+6
    prices = prices_guess
    offersprev = offers_guess
    prices_store = zeros(J,it_tols[3]+2)
    prices_store[:,2] = prices
    while (maxdiff > it_tols[4]) & (iter < it_tols[3] + 1)
        iter += 1
        # Find wage equilibrium given current guess
        offers, Z, matches = wage_eq(para, Kl, rents, hhs, hl, firms, fl, offersprev, prices, order, it_tols)
        wages =
        # Calculate marginal costs
        prices_new = (1/α).*(Z./Kl[firms[:]])^((1/α) - 1).*wages
        prices_store[:,iter+3] = prices_new
        # calculate maximum difference between prices
        maxdiff = findmax(abs.(prices - prices_new))[1]

        # Check for convergence
        if maxdiff < it_tols[4]
            println("Price equilibrium achieved after $iter iterations")
        else
            prices = 0.4*(prices_new + prices) +
            0.2*(prices_store[:,iter + 1] + prices_store[:,iter])
        end # convergence check
    end # price convergence loop
    return prices, offers, matches
end # price_eq
