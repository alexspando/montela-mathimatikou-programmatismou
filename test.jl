###############################################################
# BENDERS DECOMPOSITION FOR POWER SYSTEM EXPANSION PLANNING  #
# With LOL (unserved energy) inside the demand constraints    #
# Ensures slave is ALWAYS feasible → ONLY optimality cuts     #
###############################################################

using JuMP
using GLPK
using CSV
using DataFrames
using Printf

println("=== LOADING INPUT DATA ===")

##############################
# Load Needs Data (Slices)   #
##############################
needs = CSV.read("needs.csv", DataFrame)
println("\nNeeds:")
println(needs)

categories = needs.Charge_category
T = Float64.(needs.duration)
minL = Float64.(needs.min_level)
maxL = Float64.(needs.max_level)

m = length(categories)

ΔD = [maxL[j] - minL[j] for j in 1:m]   # width of each demand slice
println("\nΔD per slice = ", ΔD)


#######################################
# Load Technologies (MC and Inv cost) #
#######################################
tech = CSV.read("technology.csv", DataFrame)
println("\nTechnologies:")
println(tech)

MC = Float64.(tech.cost)
I  = Float64.(tech.initial_investment)
names = tech.technology
n = length(names)


#######################################
# Parameter: Value of Lost Load (VOLL)
#######################################
V = 1000.0f0   # υψηλή τιμή για να μηδενίζει το lol στο optimum


#############################################################
# MASTER PROBLEM (Investment problem with θ & Benders cuts)
#############################################################
master = Model(GLPK.Optimizer)

@variable(master, x[1:n] >= 0)   # investment in MW
@variable(master, θ >= 0)

@objective(master, Min, sum(I[i] * x[i] for i in 1:n) + θ)

cuts = Vector{ConstraintRef}()   # store Benders cuts here

println("\n=== STARTING BENDERS ITERATIONS ===")

max_iters = 50
UB = +Inf
LB = -Inf

results = DataFrame(
    iter = Int[],
    gas = Float64[],
    nuclear = Float64[],
    coal = Float64[],
    oil = Float64[],
    Q = Float64[],
    LB = Float64[],
    UB = Float64[],
    gap = Float64[],
    cut_type = String[]
)

##############################
# MAIN LOOP
##############################
for k in 1:max_iters
    println("\n---------------------------")
    println("Iteration $k")
    println("---------------------------")

    #########################
    #Solve MASTER
    #########################
    optimize!(master)
    if termination_status(master) != MOI.OPTIMAL
        println("Master not optimal. Stopping.")
        break
    end

    xk = value.(x)
    θk = value(θ)
    investment_cost = sum(I[i] * xk[i] for i in 1:n)

    println("[Master] x = ", xk)
    println("[Master] Investment cost = ", investment_cost)
    println("[Master] θ = ", θk)

    global LB
    LB = objective_value(master)


    #########################
    # 2. Solve SLAVE
    #########################
    sub = Model(GLPK.Optimizer)

    @variable(sub, p[1:n,1:m] >= 0)
    @variable(sub, lol[1:m] >= 0)

    @objective(sub, Min,
        sum(MC[i] * T[j] * p[i,j] for i in 1:n, j in 1:m) +
        sum(V     * T[j] * lol[j] for j in 1:m)
    )

    # Demand with LOL
    @constraint(sub, demand[j=1:m],
        sum(p[i,j] for i in 1:n) + lol[j] == ΔD[j]
    )

    # Capacity limits
    @constraint(sub, cap[i=1:n],
        sum(p[i,j] for j in 1:m) <= xk[i]
    )

    optimize!(sub)

    if termination_status(sub) != MOI.OPTIMAL
        println("Slave NOT optimal → THIS SHOULD NOT HAPPEN since LOL makes it feasible.")
        break
    end

    Qx = objective_value(sub)
    println("[Slave] Q(x) = ", Qx)

    global UB
    UB = min(UB, investment_cost + Qx)

    # Dual multipliers
    λ = dual.(demand)
    ρ = dual.(cap)

    println("Dual λ = ", λ)
    println("Dual ρ = ", ρ)

    gap = UB - LB
    println("Gap = ", gap)

    # Store results
    push!(results, (
        k,
        xk[1], xk[2], xk[3], xk[4],
        Qx, LB, UB, gap,
        "optimality"
    ))

    #########################
    # Check convergence
    #########################
    if abs(gap) < 1e-3
        println("\n>>> CONVERGED <<<")
        break
    end

    #########################
    # Add NEW BENDERS CUT
    #########################
    cut = @constraint(master,
        θ >= sum(λ[j] * ΔD[j] for j in 1:m) +
             sum(ρ[i] * x[i]  for i in 1:n)
    )
    push!(cuts, cut)

    println("Added optimality cut.")
end

println("\n===== FINAL SOLUTION =====")
println("x* = ", value.(x))
println("θ* = ", value(θ))
println("\nResults Table:")
println(results)

open("benders_results.csv","w") do io
    CSV.write(io, results)
end
println("\nResults saved to benders_results.csv.")

