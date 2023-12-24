import math

import matplotlib.pyplot as plt
import numpy as np

protonFileName = '../data/50000/0-10/Be_0-10.dat'


def calculate_rapidity(p0, pz):
    return 0.5 * math.log((p0 + pz) / (p0 - pz))


# 解析快度区间字符串
rapidity_range_strings = [
    "-0.1<y<0.0", "-0.2<y<-0.1", "-0.3<y<-0.2", "-0.4<y<-0.3", "-0.5<y<-0.4",
    "-0.6<y<-0.5", "-0.7<y<-0.6", "-0.8<y<-0.7", "-0.9<y<-0.8", "-1.0<y<-0.9"
]

rapidity_ranges = {}
for range_string in rapidity_range_strings:
    parts = range_string.split('<y<')
    min_y, max_y = map(float, parts)
    rapidity_ranges[range_string] = (min_y, max_y)

# 初始化存储结构
ptBins = 10
d_pt = 0.2
protonPtsByRapidity = {label: np.zeros(ptBins) for label in rapidity_ranges}
protonCountsByRapidity = {label: 0 for label in rapidity_ranges}
total_events = 0

# 读取并处理数据文件
with open(protonFileName, 'r') as file:
    lines = file.readlines()
    i = 0
    while i < len(lines):
        line = lines[i].strip()

        if line.startswith('Number of events'):
            num_events = int(line.split(': ')[1])
            i += 2  # 跳过标题行

            for j in range(num_events):
                data = lines[i + j].strip().split()
                px, py, pz, p0 = map(float, data[4:8])
                pt = math.sqrt(px ** 2 + py ** 2)
                rapidity = calculate_rapidity(p0, pz)

                for label, (min_y, max_y) in rapidity_ranges.items():
                    if min_y < rapidity < max_y:
                        npt = int(pt / d_pt)
                        if npt < ptBins:
                            protonPtsByRapidity[label][npt] += 1 / num_events ** 8
                        protonCountsByRapidity[label] += 1
                        break

            total_events += num_events
            i += num_events
        else:
            i += 1

# 绘图
plt.figure(figsize=(12, 8))
for label, pts in protonPtsByRapidity.items():
    normalized_pts = [pt / (2 * math.pi * (i * d_pt + d_pt / 2) * d_pt * total_events) for i, pt in enumerate(pts)]
    plt.plot([(i * d_pt + d_pt / 2) for i in range(ptBins)], normalized_pts, label=f'{label}')
plt.title('Normalized Proton $p_T$ Distribution by Rapidity')
plt.xlabel('$p_T$ (GeV/c)')
plt.ylabel('Normalized Yield')
plt.legend()
plt.grid(True)
plt.show()
