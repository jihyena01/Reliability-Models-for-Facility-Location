
function RPMP(I, J, h, d, NF, F, u, q, P, alpha)

    # the Reliabiltiy P-Median Problem
    m = Model(CPLEX.Optimizer)

    @variable(m, X[j in J], Bin)
    @variable(m, Y[i in I, j in J, r in (0:P-1)], Bin)
    w1 = sum(h[i] * d[i,j] * Y[i,j,0] for i in I, j in J)
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

function relaxed_RPMP(I, J, h, d, NF, F, u, q, P, alpha, lambda)

    # the Relaxed P-Median Problem
    m = Model(CPLEX.Optimizer)

    @variable(m, X[j in J], Bin)
    @variable(m, Y[i in I, j in J, r in 0:(P-1)], Bin)

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
    relaxed_term = sum(lambda[i,(r+1)] * (1-sum(Y[i,j,r] for j in J) -(r > 0 ? sum(Y[i,j,s] for j in NF, s in 0:(r-1)) : 0)) for i in I, r in 0:(P-1))
    @objective(m, Min, sum(psi[i,j,(r+1)]*Y[i,j,r] for i in I, j in J, r in 0:(P-1)) + relaxed_term)

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


function RFLP(I, J, h, d, f, NF, F, u, q, alpha)
    # P is not defined! 
    M = length(J)
    # the Reliabiltiy Fixed-Charge Location Problem
    m = Model(CPLEX.Optimizer)

    @variable(m, X[j in J], Bin)
    @variable(m, Y[i in I, j in J, r in (0:M-1)], Bin)
    w1 = sum(f[j]*X[j] for j in J) + sum(h[i] * d[i,j] * Y[i,j,0] for i in I, j in J)
    w2 = sum(h[i] * (sum(d[i,j]* q^r * Y[i,j,r] for j in NF, r in 0:(M-1)) + sum(d[i,j] * q^r * (1-q) * Y[i,j,r] for j in F, r in 0:(M-1)) ) for i in I)
    @objective(m, Min, sum(alpha * w1 + (1-alpha) * w2))

    @constraint(m, c1[i in I, r in 0:(M-1)], sum(Y[i,j,r] for j in J) + sum(Y[i,j,s] for j in NF, 
    s in 0:(r-1)) == 1)
    @constraint(m, c2[i in I, j in J, r in 0:(M-1)],Y[i,j,r] <= X[j])
 
    @constraint(m, c4[i in I, j in J], sum(Y[i,j,r] for r in 0:(M-1)) <=1)
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


function relaxed_RFLP(I, J, h, d, NF, F, u, q, alpha, lambda)
    # P is not defined! 
    M = length(J)

    # the Relaxed Fixed-charge Location Problem(RFLP)
    m = Model(CPLEX.Optimizer)

    @variable(m, X[j in J], Bin)
    @variable(m, Y[i in I, j in J, r in 0:(P-1)], Bin)

    psi = Array{Float64, 3}(undef, length(I), length(J), M) # r인덱스 1부터 시작
    for i in I
        for j in J
            for r in 0:(M-1)

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
    psi_tilde = Array{Float64, 3}(undef, length(I), length(J), M) # r인덱스 1부터 시작
    for i in I
        for j in J
            for r in 1:M
                if j in F
                    psi_tilde[i,j,r] = psi[i,j,r] - lambda[i,r]
                else # j in NF
                    psi_tilde[i,j,r] = psi[i,j,r] - (sum(lambda[i,s] for s in r : M))
                end
            end
        end
    end

    # obj_term = sum(psi_tilde[i,j,(r+1)] * Y[i,j,r] for i in I, j in J, r in 0:(M-1)) + sum(lambda[i,(r+1)] for i in I, r in 0:(M-1))
    
    w1 = sum(f[j] * X[j] for j in J) + sum(h[i] * d[i,j] * Y[i,j,0] for i in I, j in J)
    w2 = sum(h[i] * (sum(d[i,j]* q^r * Y[i,j,(r)] for j in NF, r in 0:(M-1)) + sum(d[i,j] * q^r * (1-q) * Y[i,j,(r)] for j in F, r in 0:(M-1)) ) for i in I)

    @objective(m, Min, alpha * w1 + (1-alpha) * w2)
    # @objective(m, Min, alpha * sum(f[j]*X[j] for j in J) + obj_term)

    @constraint(m, c2[i in I, j in J, r in 0:(M-1)],Y[i,j,r] <= X[j])
    @constraint(m, c3[i in I, j in J], sum(Y[i,j,r] for r in 0:(M-1)) <=1)
    @constraint(m, c4, X[u] .==1)

    optimize!(m)

    if termination_status(m) == MOI.OPTIMAL
        println("The optimal value is: ", objective_value(m))
        println("The optimal solution is: ", value.(X))
    else
        println("No solution")
    end

    return value.(X), value.(Y), objective_value(m)
end


function RFLP(f, I, J, h, d, NF, F, u, q, P, alpha)

    # the Reliabiltiy P-Median Problem
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



function PMP(I, J, h, d, P)

    # the Basic P-Median Problem
    m = Model(CPLEX.Optimizer)

    @variable(m, X[j in J], Bin)  # Decision variable for facility location
    @variable(m, Y[i in I, j in J], Bin)  # Decision variable for assignment of customers to facilities

    # Objective: Minimize the total weighted distance between facilities and customers
    @objective(m, Min, sum(h[i] * d[i,j] * Y[i,j] for i in I, j in J))

    # Constraint: Each customer is assigned to exactly one facility
    @constraint(m, c1[i in I], sum(Y[i,j] for j in J) == 1)

    # Constraint: A customer can be assigned to a facility only if the facility is open
    @constraint(m, c2[i in I, j in J], Y[i,j] <= X[j])

    # Constraint: Exactly P facilities are open
    @constraint(m, c3, sum(X[j] for j in J) == P)

    optimize!(m)

    if termination_status(m) == MOI.OPTIMAL
        println("The optimal value is: ", objective_value(m))
        println("The optimal solution is: ", value.(X))
    else
        println("No solution")
    end

    return value.(X), value.(Y), objective_value(m)
end