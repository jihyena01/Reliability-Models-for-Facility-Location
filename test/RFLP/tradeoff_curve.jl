using JuMP, CPLEX
using Random, StatsBase, Plots
import Random
using CSV, DataFrames

include("../../src/formulation.jl")
include("../../src/heuristic.jl")
include("../../data/customed_customers.jl")
include("../../data/hospital_customers.jl")
include("../../src/Lagrangian_relaxation.jl")

seed = 123
Random.seed!(seed)

## Data loading
# I, J, h, d, f, NF, F, u, q, P, alpha, lambda, _, _ = data_setting()
I, J, h, d, f, NF, F, u, q, P, alpha, lambda = hospital_data_setting()

#### trade_off_curve ---------------------------------------------------
alpha_values = [0.0, 1.0]
alpha_X_1, alpha_Y_1, obj_1 = RFLP(I, J, h, d, f, NF, F, u, q, alpha_values[1])
alpha_X_2, alpha_Y_2, obj_2 = RFLP(I, J, h, d, f, NF, F, u, q, alpha_values[2])

x_solutions = [alpha_X_1, alpha_X_2]
y_solutions = [alpha_Y_1, alpha_Y_2]
objective_values = [obj_1, obj_2]

# 절충 곡선을 DataFrame으로 관리
trade_off_curve = DataFrame(alpha = alpha_values, x_solution = x_solutions, y_solution = y_solutions, objective = objective_values)

# 절충 곡선 업데이트
solve_RPMP_and_update_curve!(trade_off_curve)

# 결과 출력
println(trade_off_curve)


# ------------------------------------------------------------------
# alpha = 0, 1일 때의 RFLP 솔루션 값이 달라지지 않음(random). // hospital은 1개 생성됨.
# 임의로 alpha 생성

alpha_values = collect(range(0.0, stop=1.0, step=0.01))

UFLP_results = []
expected_failure_cost = []

for alpha in alpha_values
    alpha_X, alpha_Y, _ = RFLP(I, J, h, d, f, NF, F, u, q, alpha)
    alpha_w1, alpha_w2 = get_multiobjective_value(alpha_X, alpha_Y, I, J, h, d, NF, F, q, P)

    push!(UFLP_results, alpha_w1)
    push!(expected_failure_cost, alpha_w2-alpha_w1)
end

plot(UFLP_results, expected_failure_cost, label=false)
scatter!(UFLP_results, expected_failure_cost, label="points")
xlabel!("UFLP Results")
ylabel!("Expected Failure Cost")
title!("Graph of Expected Failure Cost against UFLP Results)")

savefig("trade_off_curve_for_hospital_data.png")
# graph for NF = 5, 9
# alpha 값이 0, 1 일 때의 값이 매우 크거나, 음수가 되므로 해당 양끝값은 제외하였음.

# ------------------------------------------------------------------


function find_adjacent_pairs(trade_off_curve)
    adjacent_pairs = []
    for i in 1:(nrow(trade_off_curve) - 1)
        pair = (trade_off_curve[i, :], trade_off_curve[i + 1, :])
        push!(adjacent_pairs, pair)
    end
    return adjacent_pairs
end

function solve_RPMP_and_update_curve!(trade_off_curve)
    while true
        adjacent_pairs = find_adjacent_pairs(trade_off_curve)
        all_explored = true

        for pair in adjacent_pairs
            x_sol_1, y_sol1, alpha1 = pair[1][:x_solution], pair[1][:y_solution], pair[1][:alpha]
            x_sol_2, y_sol2, alpha2 = pair[2][:x_solution], pair[2][:y_solution], pair[2][:alpha]
            println(alpha1)
            println(alpha2)
            w1_1, w2_1 = get_multiobjective_value(x_sol_1, y_sol1, I, J, h, d, NF, F, q, P)
            w1_2, w2_2 = get_multiobjective_value(x_sol_2, y_sol2, I, J, h, d, NF, F, q, P)
            
            new_alpha = -(w2_1 - w2_2) / ((w1_1 - w1_2) - w2_1 + w2_2)
            
            if new_alpha > alpha1 && new_alpha < alpha2
                new_x_sol, new_y_sol, new_obj = RPMP(I, J, h, d, NF, F, u, q, P, new_alpha)
                
                if !any(row -> row[:alpha] == new_alpha, eachrow(trade_off_curve))
                    new_row = DataFrame(alpha = [new_alpha], x_solution = [new_x_sol], y_solution = [new_y_sol], objective = [new_obj])
                    append!(trade_off_curve, new_row)
                    sort!(trade_off_curve, :alpha)
                    all_explored = false
                end
            end
        end

        if all_explored
            break
        end
    end
end

function get_multiobjective_value(X, Y, I, J, h, d, NF, F, q, P)
    w1 = sum(f[j] * X[j] for j in J) + sum(h[i] * d[i,j] * Y[i,j,0] for i in I, j in J)
    w2 = sum(h[i] * (sum(d[i,j]* q^r * Y[i,j,(r)] for j in NF, r in 0:(P-1)) + sum(d[i,j] * q^r * (1-q) * Y[i,j,(r)] for j in F, r in 0:(P-1)) ) for i in I)
    return w1, w2
end


