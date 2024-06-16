using JuMP, CPLEX
using Random, StatsBase
import Random
using CSV, DataFrames

include("../../src/formulation.jl")
include("../../src/heuristic.jl")
include("../../data/hospital_customers.jl")
include("../../src/Lagrangian_relaxation.jl")

seed = 123
Random.seed!(seed)

# ------------------------------------------------------------------
# 데이터 로드
# P 개수 100개로 설정
I, J, h, d, f, NF, F, u, q, P, alpha, lambda = hospital_data_setting()
alpha = 0
# ------------------------------------------------------------------
# RPMP와 Lagrangian relaxation 알고리즘을 이용한 RPMP의 결과 비교

origin_X, origin_Y, origin_val = RFLP(I, J, h, d, f, NF, F, u, q, alpha)
X, Y, obj_val = relaxed_RFLP(I, J, h, d, NF, F, u, q, alpha, lambda)
UB, LB, heuristic_val, heuristic_X, heuristic_Y, n = Lagrangian_relaxation_RFLP(I, J, h, d, f, NF, F, u, q, P, alpha)  # heuristic_Y(r index) 1부터 시작


# ------------------------------------------------------------------
# 실험 결과 출력

epsilon = 10^-6
if abs(UB - origin_val) < epsilon
    println("The optimal solution is found by heuristic!")
else
    println("The optimal solution is not found!")
end

println("--------------------------------------------")
println("The upeer bound is: ", UB)
println("The lower bound is: ", LB)
println("The iteration number is: ", n)
println("Relaxed RFLP objective value : ", obj_val)
println("Original RFLP objective value at first: ", origin_val)


# ------------------------------------------------------------------
# 실험 설계
## a modification 부분 check! level-r 을 0~4까지로 수정. but, in my case, 시설의 개수가 많지 않고 nonfailure 시설에 배치될 확률이 높음
# 왜냐면 failure 확률 q를 높게 설정했기 때문. (이건 해당 reliability의 성능을 확인해보기 위하여 이렇게 설정함)
# 따라서 기존 방식을 그대로 사용!

# P, 즉 open facilities 시설 개수가 많을수록 convergence problem에 직면했다고 언급되어 있음.(두 개의 case 모두)


