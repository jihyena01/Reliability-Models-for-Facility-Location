using Random, StatsBase, Distances
import Random

function data_setting()
    seed = 123
    Random.seed!(seed)

    ## Data setting 
    # number of customers & facilities
    customers = 30
    facilities = 40
    nonfailable_faciltiy = 10

    # demand & cost
    demand_low = 10
    demand_high = 20

    cost_low = 20
    cost_high = 100
    penalty_cost = 150 # 임의로 지정함.

    fixed_cost_low = 100
    fixed_cost_high = 500

    # load the data
    I = collect(1:customers) # customers
    J = collect(1:facilities) # facilities

    h = Random.rand(demand_low:demand_high, length(I)) # demand of customer I
    # d = Random.rand(cost_low:cost_high, length(I), length(J)) # cost of customer I to facility J

    # 위치 데이터 생성
    customer_positions = rand(2, customers) # 고객 위치: 2차원 좌표, 각 고객마다 랜덤
    facility_positions = rand(2, facilities) # 시설 위치: 2차원 좌표, 각 시설마다 랜덤

    # 거리 행렬 계산
    d = zeros(customers, facilities)
    for i in 1:customers
        for j in 1:facilities
            d[i,j] = euclidean(customer_positions[:,i], facility_positions[:,j]) # 유클리드 거리 계산
        end
    end

    f = Random.rand(fixed_cost_low:fixed_cost_high, length(J)) # fixed cost of facility J
    theta = fill(penalty_cost, length(I)) # penalty cost
    NF = sample(J, nonfailable_faciltiy, replace=false)

    F = setdiff(J,NF)
    u = rand(NF, 1)
    println("emergency facility is ", u)
    # penalty cost 
    for i in I 
        d[i, u] .= theta[i]
        f[u] .= 0
    end


    q = 0.5 # failure probability
    alpha = 0.1
    P = facilities # P+1 = 5, 전체 facility의 수
    lambda = Random.rand(0:1, length(I), P)

    return I, J, h, d, f, NF, F, u, q, P, alpha, lambda, customer_positions, facility_positions
end