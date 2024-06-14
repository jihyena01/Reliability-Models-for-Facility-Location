using Distances
using DataFrames, CSV

# 위도 경도 거리 계산 함수
function haversine(lat1, lon1, lat2, lon2)
    radius = 6371.0  # 지구 반지름 (킬로미터)
    dlat = deg2rad(lat2 - lat1)
    dlon = deg2rad(lon2 - lon1)
    a = sin(dlat/2)^2 + cos(deg2rad(lat1)) * cos(deg2rad(lat2)) * sin(dlon/2)^2
    c = 2 * atan(sqrt(a), sqrt(1-a))
    distance = radius * c
    return distance
end


# 병원 위치 (시설)
facilities = DataFrame(
    name = ["Daejeon Veterans Hospital", "Daejeon Korean Hospital", "Daecheong Hospital", 
            "Konyang University Hospital", "Eulji University Hospital Daejeon", "Yuseong Sun Hospital",
            "Daejeon St. Mary's Hospital", "Daejeon Sun Hospital", "Chungnam National University Hospital"],
    lat = [36.4284, 36.3320, 36.3017, 36.2881, 36.3542, 36.3676, 36.3226, 36.3245, 36.3128],
    lon = [127.4277, 127.4455, 127.3731, 127.3441, 127.3784, 127.3374, 127.4188, 127.4216, 127.4075]
)

# 고객 위치 랜덤 생성 (대전 내)
customers = DataFrame(name = String[], lat = Float64[], lon = Float64[])
Random.seed!(42)
for i in 1:30
    lat = rand(36.25:0.0001:36.45)
    lon = rand(127.3:0.0001:127.5)
    push!(customers, ("Customer$i", lat, lon))
end

# 거리 행렬 초기화
distance_matrix = Array{Float64}(undef, nrow(customers), nrow(facilities))

# 거리 계산
for i in 1:nrow(customers)
    for j in 1:nrow(facilities)
        customer_location = (customers.lat[i], customers.lon[i])
        facility_location = (facilities.lat[j], facilities.lon[j])
        distance_matrix[i, j] = haversine(customer_location[1], customer_location[2], facility_location[1], facility_location[2])
    end
end

# 결과 출력
df = DataFrame(distance_matrix, :auto)
println(df)

# CSV 파일로 저장 
CSV.write("./data/customers.csv", customers)
CSV.write("./data/facilities.csv", facilities)
CSV.write("./data/distance_matrix.csv", df)