using Pkg
using LinearAlgebra
using Plots
using Distributions
using Parameters, Random, Combinatorics

@with_kw struct parameters
    I::Int64 = 5000     # Number of households
    J::Int64 = 100      # Number of firms
    L::Int64 = 20       # Number of locations
    σ1::Float64 = 0.5   # Utility parameter
    σ2::Float64 = 0.5   # Utility parameter
    ψ::Float64 = 0.25   # Work distaste parameter
    τ1::Float64 = 0     # Fixed cost of shopping
    τ2::Float64 = 0.5   # Variable cost of shopping
    γ1::Float64 = 0     # Fixed cost of commuting
    γ2::Float64 = 0.5   # Variable cost of commuting
    κ::Float64 = 5      # Non-sunk firm fixed cost
    α::Float64 = 0.7    # Production function DRTS parameter
    φ::Float64 = 0.1    # Wage adjustment increment
    u0::Float64 = 0     # Outside option utility
    λ::Float64 = 0.5    # Marginal cost of land
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

function γ(para)
    @unpack L, γ1, γ2 = para
    locs = collect(1:L)
    commcosts = γ1 .+ γ2.*abs.(locs .- transpose(locs))
    return commcosts
end # γ

function wage_eq(para, Kl, rents, hhs, hl, firms, fl, offers_guess, prices, order, Θ, it_tols)
    @unpack I, J, L, σ1, σ2, ψ, γ1, γ2, κ, α = para
    # Calculate firm demands
    Z = σ1*σ2.*exp.(-σ1.*(transpose(prices) .+ τ(para,firms)))*transpose(hl)

    # initialize while loop
    iter = 0
    maxdiff = 1e+6
    offers = offers_guess

    # while loop to obtain optimal firm offers (wage equilibrium)
    while (maxdiff > it_tols[2]) & (iter < it_tols[1] + 1)
        iter += 1
        # Compile matrix of potential ω's for each household
        ω_ij = (offers[:,:,1] - γ(para)[hhs,firms]).*offers[:,:,2] - ψ.*offers[:,:,2].^2
        # If greater than 0, hh's apply to job (i,j)
        applications = ω_ij .> 0
        # Firms choose which hh's to hire
        Ωhj = zeros(J,I)
        for j in 1:J
            # Find households that applied to j
            Ωaj = zeros(0)
            for i in 1:I
                if application[i,j] == true
                    append!(Ωaj, Int(i))
                end # appending if statement
            end # i-loop
            # Find required labor
            Nj = (Z[j]/Kl[firms[j]])^(1/α)
            # Find all possible combinations of hires
            all_plans = vcat([collect(combinations(1:size(Ωaj)[2],i)) for i=1:size(Ωaj)[2]]...)
            # Find plans that produce at least Nj
            Nj_costs = zeros(0)
            for x in 1:size(all_plans)[2]
                # If output of plan x >= Nj, append cost to Nj_costs
                if sum(Θ[all_plans[x]].*offers[all_plans[x],j,2]) >= N_j
                    append!(Nj_costs, sum(offers[all_plans[x],j,1].*offers[all_plans[x],j,2]))
                end
            end # x-loop
            # Search for minimum cost plan among Nj_plans
            if size(Nj_costs)[1] > 0
                plan_j = findmin(Nj_costs)[1]
                # Make offers
                Ωhj[j,plan_j] .= 1
            end # size check
        end # j-loop

        # Households accept best offer
        matches = zeros(J,I)
        hires = transpose(ω_ij) .* Ωhj
        for i in 1:I
            matches[findmax(hires[:,i]),i] = 1
        end # i-loop

        """ STOPPING CRITERIA"""
        if
            println("Wage Equilibrium achieved after $iter iterations")
        # If convergence is not achieved, update offers
        else
            for j in 1:J
                for i in 1:I
                    if (matches[j,i] = 1) & (Ω[j,i] = 1)
                        
                    else

                    end # End updating if statement
                end # i-loop
            end # j-loop
        end # convergence check
    end # convergence loop
    return offers, Z, matches
end # wage_eq


"""
MATRICES:
rents[1,L] -- rents for each location l
hhs[1,I] -- locations of each household (0 = outside)
firms[1,J] -- locations of each firm (0 = not operating)
offers[I,J,2] -- offers by firm j to household i in form (w,n)
prices[1,J] -- prices of goods by each firm j
hl[1,L+1] -- distribution of hh's (L+1 = outside)
fl[1,L+1] -- distribution of firms (L+1 = not operating)
Θ[I,J] -- Individual household productivity multipliers θ_i
matches[J,I] -- binary matrix matching households to firms; ie if matches[j,i] = 1, then i works for j
"""

it_tols = [200,1e-5,
200,1e-5,
200,1e-5,
200,1e-5]
