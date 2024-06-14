using Random, StatsBase
import Random

function hospital_data_setting()
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

    penalty_cost = 15 # 임의로 지정함.

    # load the data
    I = collect(1:customers) # customers
    J = collect(1:facilities) # facilities

    h = Random.rand(demand_low:demand_high, length(I)) # demand of customer I
    d = Matrix(distance_df) # distance between customer I and facility J
    theta = fill(penalty_cost, length(I)) # penalty cost

    NF = sample(J, nonfailable_faciltiy, replace=false)

    F = setdiff(J,NF)
    u = rand(NF, 1)
    # println("emergency facility is ", u)
    # penalty cost 
    for i in I 
    d[i, u] .= theta[i]
    end

    q = 0.5 # failure probability
    alpha = 0.1
    P = facilities # P+1 = 5, 전체 facility의 수
    lambda = Random.rand(0:1, length(I), P)

    return I, J, h, d, NF, F, u, q, P, alpha, lambda
end
