using Random, StatsBase
import Random

function hospital_data_setting()
    seed = 123
    Random.seed!(seed)

    ## Data Loading
    customers_df = CSV.read("data/customers.csv", DataFrame)
    facilities_df = CSV.read("data/facilities.csv", DataFrame)
    distance_df = CSV.read("data/distance_matrix.csv", DataFrame)

    # number of customers & facilities
    customers = nrow(customers_df)
    facilities = nrow(facilities_df)
    nonfailable_faciltiy = 8

    # demand & cost
    demand_low = 10
    demand_high = 20
    fixed_cost_low = 100
    fixed_cost_high = 500

    penalty_cost = 15 # 임의로 지정함.

    # load the data
    I = collect(1:customers) # customers
    J = collect(1:facilities) # facilities

    h = Random.rand(demand_low:demand_high, length(I)) # demand of customer I
    d = Matrix(distance_df) # distance between customer I and facility J
    theta = fill(penalty_cost, length(I)) # penalty cost
    f = Random.rand(fixed_cost_low:fixed_cost_high, length(J)) # fixed cost of facility J
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

    return I, J, h, d, f, NF, F, u, q, P, alpha, lambda
end
