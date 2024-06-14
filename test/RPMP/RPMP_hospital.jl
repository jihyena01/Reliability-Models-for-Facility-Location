using JuMP, CPLEX
using Random, StatsBase
import Random
using CSV, DataFrames

include("../../src/formulation.jl")
include("../../src/heuristic.jl")
include("../../src/Lagrangian_relaxation.jl")
include("../../data/hospital_customers.jl")

seed = 123
Random.seed!(seed)

# ------------------------------------------------------------------
# 데이터 로드
I, J, h, d, f, NF, F, u, q, P, alpha, lambda, _, _ = data_setting()




# ------------------------------------------------------------------
# RPMP와 Lagrangian relaxation 알고리즘을 이용한 RPMP의 결과 비교

X, Y, obj_val = relaxed_RPMP(I, J, h, d, NF, F, u, q, P, alpha, lambda)
origin_X, origin_Y, origin_val = RPMP(I, J, h, d, NF, F, u, q, P, alpha)
UB, LB, heuristic_val, heuristic_X, heuristic_Y, n = Lagrangian_relaxation(I, J, h, d, NF, F, u, q, P, alpha)  # heuristic_Y(r index) 1부터 시작




# ------------------------------------------------------------------
# 실험 결과 출력

epsilon = 10^-6
if heuristic_val - origin_val < epsilon
    println("The optimal solution is found by heuristic!")
else
    println("The optimal solution is not found!")
end


println("The upeer bound is: ", UB)
println("The lower bound is: ", LB)
println("The iteration number is: ", n)
println("Relaxed RPMP objective value : ", obj_val)
println("Original RPMP objective value at first: ", origin_val)


