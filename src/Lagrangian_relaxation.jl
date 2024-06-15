
using JuMP, CPLEX
using Random, StatsBase
import Random
using CSV, DataFrames

function Lagrangian_relaxation_RPMP(I, J, h, d, NF, F, u, q, P, alpha)
    #### Lagrangian Relaxation Algorithm ---------------------------------------------------

    ## Step 0, Initial parameter setting ==================================================

    d_var = sum(d[i,j] for i in I, j in J) / (length(I)*length(J))
    tolerance = 0.001
    beta = 2
    halving_num = 30
    min_beta = 10^(-8)
    n = 0 # iteration number
    n_max = 50 # 1200

    # Initial Multipliers
    lambda = Array{Float64, 2}(undef, length(I), P)
    for i in I
        for r in 0:(P-1)
            lambda[i,(r+1)] = h[i]*d_var / 10^(r+2)
        end
    end

    # initial upper bound & lower bound 
    UB = 10^8
    heuristic_val = -10^8
    heuristic_X, heuristic_Y = zeros(length(I), length(J), P), zeros(length(I), length(J), P)
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
            heuristic_X, heuristic_Y, heuristic_val = heuristic_RPMP(LR_X, LR_Y, I, J, P, d, q, NF, alpha)
            println("heuristic_val: ", heuristic_val)
            if heuristic_val > LB
                LB = heuristic_val
            end

            if heuristic_val < UB
                UB = heuristic_val
            end

            if heuristic_val == origin_val
                println("stop by heuristic algorithm, found optimal solution!")
                println("iteratino number: ", n)
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
        cal = n % halving_num
        if cal == 0
            beta = beta/2
        end

        if beta < min_beta
            println("stop by beta condition")
            exit_loop = true
            break #종료
        end
    end
    return UB, LB, heuristic_val, heuristic_X, heuristic_Y, n
end


function Lagrangian_relaxation_RFLP(I, J, h, d, f, NF, F, u, q, P, alpha)
    #### Lagrangian Relaxation Algorithm ---------------------------------------------------

    ## Step 0, Initial parameter setting ==================================================

    d_var = sum(d[i,j] for i in I, j in J) / (length(I)*length(J))
    tolerance = 0.001
    beta = 2
    halving_num = 30
    min_beta = 10^(-8)
    n = 0 # iteration number
    n_max = 50

    # Initial Multipliers
    lambda = Array{Float64, 2}(undef, length(I), P)
    for i in I
        for r in 0:(P-1)
            lambda[i,(r+1)] = h[i]*d_var / 10^(r+2)
        end
    end

    # initial upper bound & lower bound 
    UB = 10^8
    heuristic_val = -10^8
    heuristic_X, heuristic_Y = zeros(length(I), length(J), P), zeros(length(I), length(J), P)
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
            heuristic_X, heuristic_Y, heuristic_val = heuristic_RFLP(LR_X, LR_Y, I, J, P, f, d, q, NF, alpha)
            println("heuristic_val: ", heuristic_val)
            if heuristic_val > LB
                LB = heuristic_val
            end

            if heuristic_val < UB
                UB = heuristic_val
            end

            if heuristic_val == origin_val
                println("stop by heuristic algorithm, found optimal solution!")
                println("iteratino number: ", n)
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
        cal = n % halving_num   
        if cal == 0
            beta = beta/2
        end

        if beta < min_beta
            println("stop by beta condition")
            exit_loop = true
            break #종료
        end
    end
    return UB, LB, heuristic_val, heuristic_X, heuristic_Y, n
end