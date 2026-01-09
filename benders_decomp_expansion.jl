###############################################################
# BENDERS DECOMPOSITION FOR POWER SYSTEM EXPANSION PLANNING
###############################################################

using JuMP
using GLPK
using CSV
using DataFrames
using Printf
import MathOptInterface as MOI

function run_benders()
    println("=== LOADING INPUT DATA ===")

    ##############################
    # Load Needs Data (Slices)   #
    ##############################
    needs = CSV.read("needs.csv", DataFrame)
    
    println("\nNeeds:")
    println(needs)

    categories = needs[!, 1]
    T          = Float64.(needs[!, 2])
    minL       = Float64.(needs[!, 3])
    maxL       = Float64.(needs[!, 4])


    m = length(categories)

    # ΔD_j = max_level - min_level
    ΔD = [maxL[j] - minL[j] for j in 1:m]
    println("\nT (hours per slice) = ", T)
    println("ΔD per slice = ", ΔD)

    #######################################
    # Load Technologies (MC and Inv cost) #
    #######################################
    tech = CSV.read("technology.csv", DataFrame)

    println("\nTechnologies:")
    println(tech)

    ## Διάβασμα με θέση στήλης (robust σε BOM/hidden chars)
    names_tech = String.(tech[!, 1])      # technology
    MC         = Float64.(tech[!, 2])     # cost
    I          = Float64.(tech[!, 3])     # initial_investment

    n = length(names_tech)

    println("\nMC = ", MC)
    println("I  = ", I)
    println("Technologies = ", names_tech)


    V = 1000.0  # €/MWh 
    #############################################################
    # MASTER PROBLEM (Investment problem with θ & Benders cuts)
    #############################################################
    master = Model(GLPK.Optimizer)

    @variable(master, x[1:n] >= 0)   # investments (MW)
    @variable(master, θ >= 0)

    @objective(master, Min, sum(I[i] * x[i] for i in 1:n) + θ)

    println("\n=== STARTING BENDERS ITERATIONS ===")

    max_iters = 50
    tol = 1e-3

    LB = -Inf
    UB = +Inf

    res = DataFrame()
    res.iter = Int[]
    for nm in names_tech
        res[!, Symbol(nm)] = Float64[]
    end
    res.Q = Float64[]
    res.LB = Float64[]
    res.UB = Float64[]
    res.gap = Float64[]
    res.cut_type = String[]

    #############################################################
    # MAIN LOOP
    #############################################################
    for k in 1:max_iters
        println("\n---------------------------")
        println("Iteration $k")
        println("---------------------------")

        #########################
        # 1) Solve MASTER
        #########################
        optimize!(master)
        mst_status = termination_status(master)
        if mst_status != MOI.OPTIMAL
            println("Master not optimal. status = $mst_status. Stopping.")
            break
        end

        xk = value.(x)
        θk = value(θ)
        invest_cost = sum(I[i] * xk[i] for i in 1:n)

        @printf("[Master] obj = %.6e\n", objective_value(master))
        @printf("[Master] θ   = %.6e\n", θk)
        @printf("[Master] Inv = %.6e\n", invest_cost)
        println("[Master] x (MW):")
        for i in 1:n
            @printf("  %-10s : %.4f\n", names_tech[i], xk[i])
        end

        LB = objective_value(master)

        #########################
        # 2) Solve SLAVE
        #########################
        sub = Model(GLPK.Optimizer)

        @variable(sub, p[1:n, 1:m] >= 0)   # dispatch (MW)
        @variable(sub, lol[1:m] >= 0)      # unserved (MW)

        @objective(sub, Min,
            sum(MC[i] * T[j] * p[i,j] for i in 1:n, j in 1:m) +
            sum(V      * T[j] * lol[j] for j in 1:m)
        )

        # Demand with LOL (always feasible)
        @constraint(sub, demand[j=1:m],
            sum(p[i,j] for i in 1:n) + lol[j] == ΔD[j]
        )

        # Capacity limits coupled with xk
        @constraint(sub, cap[i=1:n],
            sum(p[i,j] for j in 1:m) <= xk[i]
        )

        optimize!(sub)
        sub_status = termination_status(sub)

        if sub_status == MOI.INFEASIBLE
            println("[Slave] INFEASIBLE (unexpected with LOL). Adding a simple feasibility cut.")
            @constraint(master, sum(x[i] for i in 1:n) >= maximum(ΔD))
            push!(res.iter, k)
            for i in 1:n
                push!(res[!, Symbol(names_tech[i])], xk[i])
            end
            push!(res.Q, NaN)
            push!(res.LB, LB)
            push!(res.UB, UB)
            push!(res.gap, UB - LB)
            push!(res.cut_type, "feasibility")
            continue
        elseif sub_status != MOI.OPTIMAL
            println("[Slave] status = $sub_status. Stopping.")
            break
        end

        Qx = objective_value(sub)
        @printf("[Slave] Q(x) = %.6e\n", Qx)

        UB = min(UB, invest_cost + Qx)
        gap = UB - LB
        @printf("[Bounds] LB = %.6e | UB = %.6e | gap = %.6e\n", LB, UB, gap)

        # Dual multipliers
        λ = dual.(demand)  # one per slice
        ρ = dual.(cap)     # one per technology

        println("Dual λ (per slice) = ", λ)
        println("Dual ρ (per tech)  = ", ρ)

        # Αποθήκευση iteration
        push!(res.iter, k)
        for i in 1:n
            push!(res[!, Symbol(names_tech[i])], xk[i])
        end
        push!(res.Q, Qx)
        push!(res.LB, LB)
        push!(res.UB, UB)
        push!(res.gap, gap)
        push!(res.cut_type, "optimality")

        #########################
        # Convergence check
        #########################
        if abs(gap) <= tol
            println("\n>>> CONVERGED (gap ≤ $(tol)) <<<")
            break
        end

        #########################
        # 3) Add Benders optimality cut
        #
        # θ ≥ Σ_j λ_j ΔD_j + Σ_i ρ_i x_i
        #
        # (Τα x_i εδώ είναι μεταβλητές του master, όχι τα xk.)
        #########################
        @constraint(master,
            θ >= sum(λ[j] * ΔD[j] for j in 1:m) +
                 sum(ρ[i] * x[i]  for i in 1:n)
        )

        println("Added optimality cut.")
    end

    println("\n===== FINAL SOLUTION =====")
    xstar = value.(x)
    println("x* (MW):")
    for i in 1:n
        @printf("  %-10s : %.4f\n", names_tech[i], xstar[i])
    end
    @printf("θ* = %.6e\n", value(θ))
    @printf("Final LB = %.6e | Final UB = %.6e | gap = %.6e\n", LB, UB, UB - LB)

    println("\nResults table:")
    println(res)

    CSV.write("benders_results.csv", res)
    println("\nSaved: benders_results.csv")
end

# Run when you include the file
run_benders()