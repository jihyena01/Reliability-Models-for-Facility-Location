using JuMP, CPLEX
using Random, StatsBase
import Random
using CSV, DataFrames

include("../../src/formulation.jl")
include("../../src/heuristic.jl")
include("../../data/customed_customers.jl")
in

seed = 123
Random.seed!(seed)

## Data loading
I, J, h, d, f, NF, F, u, q, P, alpha, lambda, _, _ = data_setting()


X, Y, obj_val = relaxed_RPMP(I, J, h, d, NF, F, u, q, P, alpha, lambda)
origin_X, origin_Y, origin_val = RPMP(I, J, h, d, NF, F, u, q, P, alpha)
println("Relaxed RPMP objective value : ", obj_val)
println("Original RPMP objective value : ", origin_val)



#### trade_off_curve

# 초기 절충 곡선 생성
alpha_values = [0.0, 1.0]
_, sol_1, obj_1 = RPMP(I, J, h, d, NF, F, u, q, P, alpha_values[1])
_, sol_2, obj_2 = RPMP(I, J, h, d, NF, F, u, q, P, alpha_values[2])

solutions = [sol_1, sol_2]
objective_values = [obj_1, obj_2]

# 절충 곡선을 DataFrame으로 관리
trade_off_curve = DataFrame(alpha = alpha_values, solution = solutions, objective = objective_values)

# 절충 곡선 업데이트
solve_RPMP_and_update_curve!(trade_off_curve)

# 결과 출력
println(trade_off_curve)


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
            sol_1, alpha1 = pair[1][:solution], pair[1][:alpha]
            sol_2, alpha2 = pair[2][:solution], pair[2][:alpha]
            println(alpha1)
            println(alpha2)
            w1_1, w2_1 = get_multiobjective_value(sol_1, I, J, h, d, NF, F, q, P)
            w1_2, w2_2 = get_multiobjective_value(sol_2, I, J, h, d, NF, F, q, P)
            
            new_alpha = -(w2_1 - w2_2) / ((w1_1 - w1_2) - w2_1 + w2_2)
            
            if new_alpha > alpha1 && new_alpha < alpha2
                _, new_sol, new_obj = RPMP(I, J, h, d, NF, F, u, q, P, new_alpha)
                
                if !any(row -> row[:alpha] == new_alpha, eachrow(trade_off_curve))
                    new_row = DataFrame(alpha = [new_alpha], solution = [new_sol], objective = [new_obj])
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

function get_multiobjective_value(Y, I, J, h, d, NF, F, q, P)
    w1 = sum(h[i] * d[i,j] * Y[i,j,0] for i in I, j in J)
    w2 = sum(h[i] * (sum(d[i,j]* q^r * Y[i,j,r] for j in NF, r in 0:(P-1)) + sum(d[i,j] * q^r * (1-q) * Y[i,j,r] for j in F, r in 0:(P-1)) ) for i in I)
    return w1, w2
end


