import matplotlib.pyplot as plt
import random
import folium
import pandas as pd
from geopy.distance import great_circle


# 랜덤 시드 고정
random.seed(42)

# 병원 위치 (시설)
facilities = [
    {"name": "Daejeon Veterans Hospital", "lat": 36.4284, "lon": 127.4277},
    {"name": "Daejeon Korean Hospital", "lat": 36.3320, "lon": 127.4455},
    {"name": "Daecheong Hospital", "lat": 36.3017, "lon": 127.3731},
    {"name": "Konyang University Hospital", "lat": 36.2881, "lon": 127.3441},
    {"name": "Eulji University Hospital Daejeon", "lat": 36.3542, "lon": 127.3784},
    {"name": "Yuseong Sun Hospital", "lat": 36.3676, "lon": 127.3374},
    {"name": "Daejeon St. Mary's Hospital", "lat": 36.3226, "lon": 127.4188},
    {"name": "Daejeon Sun Hospital", "lat": 36.3245, "lon": 127.4216},
    {"name": "Chungnam National University Hospital", "lat": 36.3128, "lon": 127.4075},
]

# CSV 파일 불러오기
customers_df = pd.read_csv('/home/comet/Documents/KAIST/Logistics-System-Optimization/Final_proj/data/customers.csv')


# 고객 데이터를 원하는 형식으로 할당
customers = []
for index, row in customers_df.iterrows():
    customer = {
        "name": row['name'],
        "lat": row['lat'],
        "lon": row['lon']
    }
    customers.append(customer)


# 플롯 초기화
plt.figure(figsize=(12, 8))
m = folium.Map(location=[36.3504119, 127.3845475], zoom_start=12)

# 시설 위치 플롯
for facility in facilities:
    plt.plot(facility["lon"], facility["lat"], 'bo', label='Hospital' if facilities.index(facility) == 0 else "")  # blue dot
    folium.Marker(
        location=[facility["lat"], facility["lon"]],
        popup=facility["name"],
        icon=folium.Icon(color='blue', icon='info-sign')
    ).add_to(m)

# 고객 위치 플롯
for customer in customers:
    plt.plot(customer["lon"], customer["lat"], 'r^', label='Customer' if customers.index(customer) == 0 else "")  # red triangle
    folium.Marker(
        location=[customer["lat"], customer["lon"]],
        popup=customer["name"],
        icon=folium.Icon(color='red', icon='info-sign')
    ).add_to(m)

# 고객을 가장 가까운 시설에 연결
for customer in customers:
    customer_location = (customer["lat"], customer["lon"])
    nearest_facility = min(facilities, key=lambda f: great_circle(customer_location, (f["lat"], f["lon"])).meters)
    facility_location = (nearest_facility["lat"], nearest_facility["lon"])
    plt.plot([customer["lon"], nearest_facility["lon"]], [customer["lat"], nearest_facility["lat"]], 'k-', alpha=0.5)
    folium.PolyLine(
        locations=[customer_location, facility_location],
        color='black'
    ).add_to(m)

# 범례 추가 (우측 하단)
plt.legend(loc='lower right')

# 플롯 설정
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.title('Hospital and Customer Locations in Daejeon')
plt.grid(True)

# 플롯 저장 및 표시
plt.savefig("Hospital_customer_plot_daejeon.png")
plt.show()

# 지도 저장
m.save("hospital_customer_map.html")