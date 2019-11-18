using Pkg
using LinearAlgebra
using Plots
using Distributions
using Parameters, Random, Combinatorics

@with_kw struct parameters
    I::Int64 = 100      # Number of households
    J::Int64 = 10     # Number of firms
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
    tcosts = τ1 .+ τ2.*abs.(locs .- firms)
    return tcosts
end #τ

function γ(para)
    @unpack L, γ1, γ2 = para
    locs = collect(1:L)
    commcosts = γ1 .+ γ2.*abs.(locs .- transpose(locs))
    return commcosts
end # γ

function wage_eq(para, Kl, hhs, hl, firms, offers_guess, prices, Θ, it_tols)
    @unpack I, J, L, σ1, σ2, ψ, γ1, γ2, κ, α, φ = para
    # Calculate firm demands
    Z = σ1*σ2.*exp.(-σ1.*(prices .+ τ(para,firms)))*hl[1:L]

    # initialize while loop
    iter = 0
    maxdiff = 1e+6
    offers = offers_guess
    # Define matricies to hold results
    Ωhj = zeros(J,I)        # Hire offers
    matches = zeros(J,I)    # Matches between i's and j's
    Ω_star = []             # i's actually hired
    matches_stor = zeros(J,I,Int(it_tols[1])+1)
    offers_stor = zeros(J,I,2,Int(it_tols[1])+1)
    offers_stor[:,:,:,1] = offers_guess
    # while loop to obtain optimal firm offers (wage equilibrium)
    while (maxdiff > it_tols[2]) & (iter < it_tols[1])
        iter += 1
        # Reset matricies to hold results
        Ωhj = zeros(J,I)        # Hire offers
        matches = zeros(J,I)    # Matches between i's and j's
        Ω_star = []             # i's actually hired
        # Compile matrix of potential ω's for each household
        ω_ij = (offers[:,:,1] - γ(para)[firms,hhs]).*offers[:,:,2] - ψ.*offers[:,:,2].^2
        # If greater than 0, hh's apply to job (i,j)
        applications = ω_ij .> 0

        # Find required labor impacts
        Nj = (Z./Kl[firms]).^(1/α)
        for j in 1:J
            # Find households that applied to j
            Ωaj = []
            for i in 1:I
                if applications[j,i] == true
                    append!(Ωaj, Int(i))
                end # appending if statement
            end # i-loop

            # Find all possible combinations of hires
            if size(Ωaj)[1] > 0
                # Take 10 most productive households for j
                all_plans = []
                if size(Ωaj)[1] > 10
                    all_plans = collect(combinations(sortperm(Θ[j,:])[(size(Θ[j,:])[1] -9):end]))
                else
                    allplans = collect(combinations(Ωaj))
                end # size of Ωaj check

                # Find plans that produce at least Nj
                Nj_costs = []
                for x in 1:size(all_plans)[1]
                    # If output of plan x >= Nj, append cost to Nj_costs
                    if sum(Θ[all_plans[x]].*log.(offers[j,all_plans[x],2])) >= Nj[j]
                        push!(Nj_costs, sum(offers[j,all_plans[x],1].*offers[j,all_plans[x],2]))
                    end
                end # x-loop
                # Search for minimum cost plan among Nj_plans
                if size(Nj_costs)[1] > 0
                    plan_j = findmin(Nj_costs)[2]
                    # Make offers
                    Ωhj[j,all_plans[plan_j]] .= 1
                end # end Nj_plans size check
            end # size check
        end # j-loop

        # Households accept best offer
        hire_offers = abs.(ω_ij .* Ωhj)
        for i in 1:I
            if findmax(hire_offers[:,i])[1] > 0
                matches[findmax(hire_offers[:,i])[2],i] = 1
            end
        end # i-loop

        # Create vector to store households matched to each j

        for j in 1:J
            push!(Ω_star, findall(i -> i == 1,matches[j,:]))
        end # j-loop

        # Update offers
        for j in 1:J
            for i in 1:I
                if (matches[j,i] == 1) & (Ωhj[j,i] == 1)
                    # Set optimal offer
                    offers[j,i,2] = exp.((Nj[j] - transpose(Θ[j,Ω_star[j]])*log.((offers[j,i,1]/Θ[j,i]).*(Θ[j,Ω_star[j]]./offers[j,Ω_star[j],1])))/sum(Θ[j,Ω_star[j]]))
                else
                    offers[j,i,1] += φ
                end # End updating if statement
            end # i-loop
        end # j-loop

        """ STOPPING CRITERIA"""
        if offers_stor[:,:,2,iter] == offers[:,:,2]
            println("Wage Equilibrium achieved after $iter iterations")
            CONT = false
            matches_stor[:,:,iter+1] = matches
            offers_stor[:,:,:,iter+1] = offers
        # If convergence is not achieved, update offers
        else
            # Store matches/offers from this iteration
            println("Iteration: $iter")
            matches_stor[:,:,iter+1] = matches
            offers_stor[:,:,:,iter+1] = offers

        end # convergence check
    end # convergence loop
    return offers_stor, Z, matches_stor, Ω_star, Ωhj
end # wage_eq

it_tols = [1000,1e-5, 200,1e-5, 200,1e-5, 200,1e-5]

offers_guess = zeros(J,I,2)
offers_guess[:,:,1] .= 1
offers_guess[:,:,2] .= 2

offers_ret, Z_ret, matches_ret, Ostar, Oh = wage_eq(para, Kl, hhs, hl, firms, offers_guess, prices, Θ, it_tols)

"""
MATRICES:
rents[L] -- rents for each location l
hhs[I] -- locations of each household (0 = outside)
firms[J] -- locations of each firm (0 = not operating)
offers[J,I,2] -- offers by firm j to household i in form (w,n)
prices[J] -- prices of goods by each firm j
hl[L+1] -- distribution of hh's (L+1 = outside)
fl[L+1] -- distribution of firms (L+1 = not operating)
Θ[J,I] -- Individual household productivity multipliers θ_i
matches[J,I] -- binary matrix matching households to firms; ie if matches[j,i] = 1, then i works for j
"""
@unpack I, J, L, σ1, σ2, ψ, γ1, γ2, κ, α, φ = para

Θ = rand(Uniform(1.5,2), J,I)
Kl = ones(L)
hhs = rand(1:L, I)
hl = zeros(L+1)
for l in 1:L
    hl[l] = count(i -> i == l, hhs)
end
firms = rand(1:L, J)
prices = ones(J)
