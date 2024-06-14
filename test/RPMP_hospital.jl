using JuMP, CPLEX
using Random, StatsBase
import Random
using CSV, DataFrames

include("../src/formulation.jl")
include("../src/heuristic.jl")

seed = 123
Random.seed!(seed)

## Data Loading
customers_df = CSV.read("/home/comet/Documents/KAIST/Logistics-System-Optimization/Final_proj/data/customers.csv", DataFrame)
facilities_df = CSV.read("/home/comet/Documents/KAIST/Logistics-System-Optimization/Final_proj/data/facilities.csv", DataFrame)
distance_df = CSV.read("/home/comet/Documents/KAIST/Logistics-System-Optimization/Final_proj/data/distance_matrix.csv", DataFrame)

# number of customers & facilities
customers = nrow(customers_df)
facilities = nrow(facilities_df)
nonfailable_faciltiy = 5

# demand & cost
demand_low = 10
demand_high = 20

penalty_cost = 10 # 임의로 지정함.

# load the data
I = collect(1:customers) # customers
J = collect(1:facilities) # facilities

h = Random.rand(demand_low:demand_high, length(I)) # demand of customer I

d = Matrix(distance_df) # distance between customer I and facility J
theta = fill(penalty_cost, length(I)) # penalty cost

NF = sample(J, nonfailable_faciltiy, replace=false)

F = setdiff(J,NF)
u = rand(NF, 1)
println("emergency facility is ", u)
# penalty cost 
for i in I 
    d[i, u] .= theta[i]
end


q = 0.5 # failure probability
alpha = 0.1
P = facilities # P+1 = 5, 전체 facility의 수
lambda = Random.rand(0:1, length(I), P)
X, Y, obj_val = relaxed_RPMP(I, J, h, d, NF, F, u, q, P, alpha, lambda)
origin_X, origin_Y, origin_val = RPMP(I, J, h, d, NF, F, u, q, P, alpha)
println("Relaxed RPMP objective value : ", obj_val)
println("Original RPMP objective value : ", origin_val)


## 기본 PMP_LR obj function 값 계산 => obj_val
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




#### Lagrangian Relaxation Algorithm ---------------------------------------------------

## Step 0, Initial parameter setting ==================================================

d_var = sum(d[i,j] for i in I, j in J) / (length(I)*length(J))
tolerance = 0.001
beta = 2
halving_num = 30
min_beta = 10^(-8)
n = 0 # iteration number
n_max = 1200

# Initial Multipliers
lambda = Array{Float64, 2}(undef, length(I), P)
for i in I
    for r in 0:(P-1)
        lambda[i,(r+1)] = h[i]*d_var / 10^(r+2)
    end
end

# initial upper bound & lower bound 
UB = 10^6
LR_X, LR_Y, LR_val = relaxed_RPMP(I, J, h, d, NF, F, u, q, P, alpha, lambda)
println("initianl LR_val: ", LR_val)
LB = LR_val
n = n +1
exit_loop = false

# Initial step size(t)
t = beta *(UB - LR_val) / (sum(1- sum(Y[i,j,r] for j in J) - (r > 0 ? sum(Y[i,j,s] for j in NF, s in 0:(r-1)) : 0) for i in I, r in 0:(P-1)))^2



check_feasibility = all(sum(LR_Y[i,j,r] for j in J) + (r > 0 ? sum(LR_Y[i,j,s] for j in NF, s in 0:(r-1)) : 0 ) == 1 for i in I, r in 0:(P-1))

if check_feasibility == true
    println("The initial solution is feasible. stop! ") # 종료
    UB = LR_val
    exit_loop = true
else
    println("The initial solution is infeasible. go to step 2") # Step 2로 이동료
    # multipliers update
    for i in I
        for r in 0:(P-1)
           
            lambda[i, (r+1)] = lambda[i,(r+1)] + t*(1- sum(LR_Y[i,j,r] for j in J) - (r > 0 ? sum(LR_Y[i,j,s] for j in NF, s in 0:(r-1)) : 0))
        end
    end
end

while exit_loop != true
    ## Step 1,===============================================================
    # Initial LR_val feasibility check
    LR_X, LR_Y, LR_val = relaxed_RPMP(I, J, h, d, NF, F, u, q, P, alpha, lambda)
    if LR_val > LB
        LB = LR_val
    end

    check_feasibility = all(sum(LR_Y[i,j,r] for j in J) + (r > 0 ? sum(LR_Y[i,j,s] for j in NF, s in 0:(r-1)) : 0 ) == 1 for i in I, r in 0:(P-1))
    n = n+1 

    if check_feasibility == true
        println("The initial solution is feasible. ") # 종료
        UB = LR_val
    else
        # println("The initial solution is infeasible. go to step 2") # Step 2로 이동
        # multipliers update
        for i in I
            for r in 0:(P-1)
                lambda[i, (r+1)] = lambda[i,(r+1)] + t*(1- sum(LR_Y[i,j,r] for j in J) - (r > 0 ? sum(LR_Y[i,j,s] for j in NF, s in 0:(r-1)) : 0))
            end
        end
        ## heuristic algorithm => find feasible solution, and update lower bound!
        heuristic_X, heuristic_Y, heuristic_val = heuristic_RPMP(LR_X, I, J, P, d, q, NF, alpha)
        println("heuristic_val: ", heuristic_val)
        if heuristic_val > LB
            LB = heuristic_val
        end

        if heuristic_val == origin_val
            println("stop by heuristic algorithm, found optimal solution!")
            exit_loop = true
            break
        end
        
        
    end


    # Step2 ================================================================
    # LR termination condition

    # LR_solution이 음수인 경우 발생 !!!
    if abs(UB - LR_val) / abs(LR_val) < tolerance
        println("stop by tolerance")
        println(LR_val)
        exit_loop = true
        break # 종료
    end

    if n > n_max
        println("stop by max_iteration")
        exit_loop = true
        break # 종료
    end

    # beta size update
    beta = n % 30
    if beta == 0
        beta = beta/2
    end

    if beta < min_beta
        println("stop by beta condition")
        exit_loop = true
        break #종료
    end
end

println("The upeer bound is: ", UB)
println("The lower bound is: ", LB)

if heuristic_Y == origin_Y
    println("The optimal solution is found!")
else
    println("The optimal solution is not found!")
end


println("Relaxed RPMP objective value : ", obj_val)
println("Original RPMP objective value at first: ", origin_val)


