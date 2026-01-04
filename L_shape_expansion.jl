###############################################################
# L-SHAPED METHOD FOR STOCHASTIC SYSTEM EXPANSION (LP)
# Output stored in CSV with decision evolution per iteration
###############################################################

using JuMP
using GLPK
using CSV
using DataFrames
using Printf

println("=== LOADING INPUT DATA ===")

#######################################
# Load Technologies (MC and Inv cost) #
#######################################
tech = CSV.read("technology.csv", DataFrame)
println("\nTechnologies:")
println(tech)

names = String.(tech.technology)
MC = Float64.(tech.cost)
I  = Float64.(tech.initial_investment)
n = length(names)

#######################################
# Parameter: Value of Lost Load (VOLL)
#######################################
V = 1000.0

#######################################
# Scenarios from Table 7 (2 scenarios)
#######################################
T = [1.0, 0.79908675799, 0.17123287671]  # durations (not stochastic)
m = length(T)

P = [0.10, 0.90]              # probabilities
Ω = length(P)

# ΔD[j, ω] = max-min (width of slice j in scenario ω)
ΔD = [
    7086.0  3919.0;  # base
    1918.0  3410.0;  # medium
    2165.0  2986.0   # peak
]  # size m x Ω

println("\nScenario probabilities P = ", P)
println("ΔD per (slice j, scenario ω):")
for ω in 1:Ω
    println("  ω=$ω  ΔD = ", ΔD[:, ω])
end

#############################################################
# MASTER PROBLEM
#############################################################
master = Model(GLPK.Optimizer)

@variable(master, x[1:n] >= 0)   # investment in MW
@variable(master, θ >= 0)

@objective(master, Min, sum(I[i] * x[i] for i in 1:n) + θ)

cuts = Vector{ConstraintRef}()

#############################################################
# SLAVE SOLVER FOR ONE SCENARIO ω
#############################################################
function solve_slave_one_scenario(xk::Vector{Float64}, ω::Int,
                                  MC::Vector{Float64}, T::Vector{Float64},
                                  V::Float64, ΔD::Matrix{Float64})
    n = length(MC)
    m = length(T)

    sub = Model(GLPK.Optimizer)
    set_silent(sub)

    @variable(sub, p[1:n, 1:m] >= 0)
    @variable(sub, lol[1:m] >= 0)

    @objective(sub, Min,
        sum(MC[i] * T[j] * p[i, j] for i in 1:n, j in 1:m) +
        sum(V     * T[j] * lol[j]  for j in 1:m)
    )

    @constraint(sub, demand[j=1:m],
        sum(p[i, j] for i in 1:n) + lol[j] == ΔD[j, ω]
    )

    @constraint(sub, cap[i=1:n],
        sum(p[i, j] for j in 1:m) <= xk[i]
    )

    optimize!(sub)
    if termination_status(sub) != MOI.OPTIMAL
        error("Slave not optimal for scenario ω=$ω (should not happen with LOL).")
    end

    Qω = objective_value(sub)
    λω = dual.(demand)
    ρω = dual.(cap)

    return Qω, λω, ρω
end

#############################################################
# RESULTS TABLE (decision evolution)
#############################################################
# We build it dynamically so the CSV has real tech names.
results = DataFrame()
results.iter = Int[]
for nm in names
    results[!, Symbol(nm)] = Float64[]   # x per tech
end
results.θ = Float64[]
results.investment_cost = Float64[]


results.EQ = Float64[]   # expected recourse
results.gap = Float64[]
results.cut_type = String[]

#############################################################
# MAIN L-SHAPED LOOP
#############################################################
println("\n=== STARTING L-SHAPED ITERATIONS ===")

max_iters = 50
UB = +Inf
LB = -Inf

for k in 1:max_iters
    println("\n---------------------------")
    println("Iteration $k")
    println("---------------------------")

    #########################
    # 1. Solve MASTER
    #########################
    optimize!(master)
    if termination_status(master) != MOI.OPTIMAL
        println("Master not optimal. Stopping.")
        break
    end

    xk = value.(x)
    θk = value(θ)
    investment_cost = sum(I[i] * xk[i] for i in 1:n)
    LB = objective_value(master)

    println("[Master] x = ", xk)
    println("[Master] Investment cost = ", investment_cost)
    println("[Master] θ = ", θk)
    println("[Master] LB = ", LB)

    #########################
    # 2. Solve ALL SLAVES => Expected recourse + aggregated cut
    #########################
    Qω_vals = zeros(Ω)
    Qbar = 0.0
    α = 0.0
    β = zeros(n)

    for ω in 1:Ω
        Qω, λω, ρω = solve_slave_one_scenario(xk, ω, MC, T, V, ΔD)

        Qω_vals[ω] = Qω
        Qbar += P[ω] * Qω

        α += P[ω] * sum(λω[j] * ΔD[j, ω] for j in 1:m)
        for i in 1:n
            β[i] += P[ω] * ρω[i]
        end
    end

    println("[Slaves] Qω = ", Qω_vals, "   E[Q(x)] = ", Qbar)

    #########################
    # 3. Update UB and gap
    #########################
    global UB
    UB = min(UB, investment_cost + Qbar)
    gap = UB - LB

    println("[Bounds] UB = ", UB)
    println("[Bounds] Gap = ", gap)

    #########################
    # 4. Store results row
    #########################
    row = Dict{Symbol, Any}()
    row[:iter] = k

    for i in 1:n
        row[Symbol(names[i])] = xk[i]
    end

    row[:θ] = θk
    row[:investment_cost] = investment_cost


    row[:EQ] = Qbar

    row[:gap] = gap
    row[:cut_type] = "optimality"

    push!(results, row)

    # Optional: write CSV each iteration (so you can watch it live)
    CSV.write("lshaped_results.csv", results)

    #########################
    # 5. Convergence check
    #########################
    if abs(gap) < 1e-3
        println("\n>>> CONVERGED <<<")
        break
    end

    #########################
    # 6. Add NEW L-SHAPED CUT
    #########################
    cut = @constraint(master,
        θ >= α + sum(β[i] * x[i] for i in 1:n)
    )
    push!(cuts, cut)
    println("Added aggregated optimality cut.")
end

println("\n===== FINAL SOLUTION =====")
println("x* = ", value.(x))
println("θ* = ", value(θ))
println("\nResults saved to lshaped_results.csv.")

#############################################################
# Display from CSV (as you asked)
#############################################################
println("\n=== DISPLAYING OUTPUT FROM CSV ===")
df = CSV.read("lshaped_results.csv", DataFrame)
println(df)
