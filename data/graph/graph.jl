using JuMP, CPLEX
using Random, StatsBase
import Random
using CSV, DataFrames
using Plots, Colors
using PlotlyBase, PlotlyJS

include("../src/formulation.jl")
include("../src/heuristic.jl")
include("customed_customers.jl")

seed = 123
Random.seed!(seed)

## Data loading
## This is case of random customers and facilities
### customer = 30, facilities = 10, Nonfailable facilites = 3

I, J, h, d, f, NF, F, u, q, P, alpha, lambda, customer_position, facility_position = data_setting()
println(customer_position)

X, Y, obj_val = RPMP(I, J, h, d, NF, F, u, q, P, alpha)
origin_X, origin_Y, origin_val = PMP(I, J, h, d, P)
println("Relaxed RPMP objective value : ", obj_val)
println("Original RPMP objective value : ", origin_val)

## Plotting ----------------------------------------------------------------
# For PMP solution

gr() 

# 범례를 활성화하고 위치를 조정합니다.
plt = plot(size=(800, 600), legend=(0.85, 0.6))

# 고객 위치 표시, 첫 번째 점에만 범례 추가
scatter!([customer_position[1, 1]], [customer_position[2, 1]], color=:blue, label="Customer")
for i in 2:size(customer_position, 2)
    scatter!([customer_position[1, i]], [customer_position[2, i]], color=:blue, label=false)
end

# 시설 위치 표시, 첫 번째 점에만 범례 추가
scatter!([facility_position[1, 1]], [facility_position[2, 1]], color=:red, label="Facility")
for j in 2:size(facility_position, 2)
    scatter!([facility_position[1, j]], [facility_position[2, j]], color=:red, label=false)
end

# 연결선 표시, 첫 번째 선에만 범례 추가
legend_added_for_lines = false
for i in 1:size(customer_position, 2)
    for j in 1:size(facility_position, 2)
        if origin_Y[i,j] > 0
            if !legend_added_for_lines
                plot!([customer_position[1, i], facility_position[1, j]], [customer_position[2, i], facility_position[2, j]], color=:black, label="Connection")
                legend_added_for_lines = true
            else
                plot!([customer_position[1, i], facility_position[1, j]], [customer_position[2, i], facility_position[2, j]], color=:black, label=false)
            end
        end
    end
end

display(plt)

savefig(plt, "PMP_graph.png")

# ------------------------------------------------------------------------
# For RPMP solution

plt = plot(size=(800, 600), legend=(0.85, 0.6))

# 고객 위치 표시, 첫 번째 점에만 범례 추가
scatter!([customer_position[1, 1]], [customer_position[2, 1]], color=:blue, label="Customer")
for i in 2:size(customer_position, 2)
    scatter!([customer_position[1, i]], [customer_position[2, i]], color=:blue, label=false)
end

# 시설 위치 표시, 첫 번째 점에만 범례 추가
scatter!([facility_position[1, 1]], [facility_position[2, 1]], color=:red, label="Facility")
for j in 2:size(facility_position, 2)
    scatter!([facility_position[1, j]], [facility_position[2, j]], color=:red, label=false)
end

# 실선과 점선 표시, 첫 번째 선에만 범례 추가
legend_added_for_lines = Dict(:solid => false, :dotted => false)
for i in 1:size(customer_position, 2)
    for j in 1:size(facility_position, 2)
        for r in 0:(P-1)
            if Y[i,j,r] > 0
                if r == 0
                    if !legend_added_for_lines[:solid]
                        plot!([customer_position[1, i], facility_position[1, j]], [customer_position[2, i], facility_position[2, j]], color=:black, label="r=0")
                        legend_added_for_lines[:solid] = true
                    else
                        plot!([customer_position[1, i], facility_position[1, j]], [customer_position[2, i], facility_position[2, j]], color=:black, label=false)
                    end
                else
                    if !legend_added_for_lines[:dotted]
                        plot!([customer_position[1, i], facility_position[1, j]], [customer_position[2, i], facility_position[2, j]], color=:lightgrey, linestyle=:dot, label="r>0")
                        legend_added_for_lines[:dotted] = true
                    else
                        plot!([customer_position[1, i], facility_position[1, j]], [customer_position[2, i], facility_position[2, j]], color=:lightgrey, linestyle=:dot, label=false)
                    end
                end
            end
        end
    end
end

display(plt)

savefig(plt, "RPMP_graph.png")