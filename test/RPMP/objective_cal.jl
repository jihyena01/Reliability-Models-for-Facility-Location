using JuMP, CPLEX
using Random, StatsBase
import Random
using CSV, DataFrames

include("../../src/formulation.jl")
include("../../src/heuristic.jl")
include("../../data/customed_customers.jl")

seed = 123
Random.seed!(seed)

## Data loading
I, J, h, d, f, NF, F, u, q, P, alpha, lambda, _, _ = data_setting()


X, Y, obj_val = relaxed_RPMP(I, J, h, d, NF, F, u, q, P, alpha, lambda)
origin_X, origin_Y, origin_val = RPMP(I, J, h, d, NF, F, u, q, P, alpha)
println("Relaxed RPMP objective value : ", obj_val)
println("Original RPMP objective value : ", origin_val)


## 기본 RPMP_LR obj function 값 계산 => obj_val
psi = Array{Float64, 3}(undef, length(I), length(J), P) # r인덱스 1부터 시작 
for i in I
    for j in J
        for r in 0:(P-1)
            if r == 0 & j in NF
                psi[i,j,(r+1)] = h[i]*d[i,j]
            elseif r == 0 & j in F
                psi[i,j,(r+1)] = alpha*h[i]*d[i,j] + (1-alpha)*h[i]*d[i,j]*(1-q)
            elseif r > 0 & j in NF
                psi[i,j,(r+1)] = (1-alpha)*h[i]*d[i,j]*q^r
            else
                psi[i,j,(r+1)] = alpha*h[i]*d[i,j]*q^r*(1-q)
            end
        end
    end
end
relaxed_term = sum(lambda[i,(r+1)] * (1- sum(Y[i,j,r] for j in J) - (r > 0 ? sum(Y[i,j,s] for j in NF, s in 0:(r-1)) : 0)) for i in I, r in 0:(P-1))
obj_val = sum(psi[i,j,(r+1)]*Y[i,j,r] for i in I, j in J, r in 0:(P-1)) + relaxed_term


## psi_tilde로 변형한 형태의 obj_function 값 계산 => obj_val_2
psi_tilde = Array{Float64, 3}(undef, length(I), length(J), P) # r인덱스 1부터 시작
for i in I
    for j in J
        for r in 1:P
            if j in F
                psi_tilde[i,j,r] = psi[i,j,r] - lambda[i,r]
            else # j in NF
                psi_tilde[i,j,r] = psi[i,j,r] - (sum(lambda[i,s] for s in r : P))
            end
        end
    end
end

obj_val_2 = sum(psi_tilde[i,j,(r+1)] * Y[i,j,r] for i in I, j in J, r in 0:(P-1)) + sum(lambda[i,(r+1)] for i in I, r in 0:(P-1))

gamma = Array{Float64, 1}(undef, length(J))
for j in J
    gamma[j] = sum( min(0, minimum(collect(psi_tilde[i,j,r] for r in 1:P))) for i in I)
end
println("gamma: ", gamma)

obj_val_3 = sum(gamma[j] * X[j] for j in J) + sum(lambda[i,(r+1)] for i in I, r in 0:(P-1))

println("obj_val: ", obj_val)
println("obj_val_2: ", obj_val_2)
println("obj_val_3: ", obj_val_3)
