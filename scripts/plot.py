import matplotlib.pyplot as plt
import numpy as np

# 定义质子的静止质量
proton_mass = 0.938  # GeV/c^2


# 读取数据
def read_data(file_name):
    with open(file_name, 'r') as file:
        lines = file.readlines()

    protons, event_count = [], 0
    for line in lines:
        if line.startswith('t'):
            event_count += 1
        else:
            parts = line.split()
            px, py, pz, p0 = map(float, parts[4:8])
            protons.append([px, py, pz, p0])
    return protons, event_count


# 计算快度
def calculate_rapidity(pz, p0):
    return 0.5 * np.log((p0 + pz) / (p0 - pz))


# 读取文件
protons, total_events = read_data('../data/100/deuteron_no_centrality.dat')

# 定义快度范围
rapidity_ranges = {
    "-1.0<y<-0.9": (-1.0, -0.9),
    "-0.9<y<-0.8": (-0.9, -0.8),
    "-0.8<y<-0.7": (-0.8, -0.7),
    "-0.7<y<-0.6": (-0.7, -0.6),
    "-0.6<y<-0.5": (-0.6, -0.5),
    "-0.5<y<-0.4": (-0.5, -0.4),
    "-0.4<y<-0.3": (-0.4, -0.3),
    "-0.3<y<-0.2": (-0.3, -0.2),
    "-0.2<y<-0.1": (-0.2, -0.1),
    "-0.1<y<0.0": (-0.1, 0.0),
    # ... 其他快度区间 ...
}

# 初始化计数器
proton_counts = {label: 0 for label in rapidity_ranges.keys()}

# 处理数据
for proton in protons:
    px, py, pz, p0 = proton
    rapidity = calculate_rapidity(pz, p0)
    for label, (min_y, max_y) in rapidity_ranges.items():
        if min_y <= rapidity < max_y:
            proton_counts[label] += 1

# 归一化计数
average_proton_counts = {label: count / total_events for label, count in proton_counts.items()}

# 输出结果
print("Average number of protons per event in each rapidity range:")
for label, avg_count in average_proton_counts.items():
    print(f"{label}: {avg_count:.2f}")

# 可视化结果
plt.figure(figsize=(10, 6))
labels = list(average_proton_counts.keys())
values = list(average_proton_counts.values())
plt.bar(labels, values)
plt.xlabel('Rapidity Range')
plt.ylabel('Average Proton Count')
plt.title('Average Proton Count per Event in Different Rapidity Ranges')
plt.xticks(rotation=45)
plt.show()
