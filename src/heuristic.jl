
function find_indices_of_ones(dense_array::JuMP.Containers.DenseAxisArray)
    indices_of_ones = Int[]
    for index in eachindex(dense_array)
        if dense_array[index] == 1.0
            push!(indices_of_ones, index[1]) 
        end
    end
    return indices_of_ones
end



function heuristic_RPMP(LR_X, I, J, P, d, q, NF, alpha)
    heuristic_X = find_indices_of_ones(LR_X)
    heuristic_Y = zeros(Int, length(I), length(J), P) # r index 1부터 시작!

    for i in I
        distance_vector = [d[i, j] for j in J]
        for j in J
            if !(j in heuristic_X)
                distance_vector[j] = 1000000
            end
        end
        for r in 1:P
            min_distance = minimum(distance_vector)
            min_index = findfirst(==(min_distance), distance_vector)
            heuristic_Y[i, min_index, r] = 1
            # println(i, min_index, r)
            distance_vector[min_index] = 1000000

            if min_index in NF
                break
            end
        end
    end
    w1 = sum(h[i] * d[i,j] * heuristic_Y[i,j,1] for i in I, j in J)
    w2 = sum(h[i] * (sum(d[i,j]* q^r * heuristic_Y[i,j,(r+1)] for j in NF, r in 0:(P-1)) + sum(d[i,j] * q^r * (1-q) * heuristic_Y[i,j,(r+1)] for j in F, r in 0:(P-1)) ) for i in I)
    heuristic_val = alpha * w1 + (1-alpha) * w2
    return LR_X, heuristic_Y, heuristic_val
end
