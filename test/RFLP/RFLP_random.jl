using JuMP, CPLEX
using Random, StatsBase
import Random
using CSV, DataFrames

include("../../src/formulation.jl")
include("../../src/heuristic.jl")
include("../../data/customed_customers.jl")

I, J, h, d, f, NF, F, u, q, P, alpha, lambda = data_setting() 



function RFLP(I, J, h, d, f, NF, F, u, q, alpha)
    # P is not defined! 

    # the Reliabiltiy Fixed-Charge Location Problem
    m = Model(CPLEX.Optimizer)

    @variable(m, X[j in J], Bin)
    @variable(m, Y[i in I, j in J, r in (0:P-1)], Bin)
    w1 = sum(f[j]*X[j] for j in J) + sum(h[i] * d[i,j] * Y[i,j,0] for i in I, j in J)
    w2 = sum(h[i] * (sum(d[i,j]* q^r * Y[i,j,r] for j in NF, r in 0:(P-1)) + sum(d[i,j] * q^r * (1-q) * Y[i,j,r] for j in F, r in 0:(P-1)) ) for i in I)
    @objective(m, Min, sum(alpha * w1 + (1-alpha) * w2))

    @constraint(m, c1[i in I, r in 0:(P-1)], sum(Y[i,j,r] for j in J) + sum(Y[i,j,s] for j in NF, 
    s in 0:(r-1)) == 1)
    @constraint(m, c2[i in I, j in J, r in 0:(P-1)],Y[i,j,r] <= X[j])
    @constraint(m, c3, sum(X[j] for j in J) == P)
    @constraint(m, c4[i in I, j in J], sum(Y[i,j,r] for r in 0:(P-1)) <=1)
    @constraint(m, c5, X[u] .==1)

    optimize!(m)

    if termination_status(m) == MOI.OPTIMAL
        println("The optimal value is: ", objective_value(m))
        println("The optimal solution is: ", value.(X))
    else
        println("No solution")
    end

    return value.(X), value.(Y), objective_value(m)
end