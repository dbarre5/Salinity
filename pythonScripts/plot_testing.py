import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.pyplot import figure
df = pd.read_csv("contour.csv")
x = df["Points:0"]*100
z = df["Points:2"]*100


#Read data for vertical lines
df = pd.read_csv("test.csv")
velList = []
vertLocList = []
horLocList = []
velList.append(df["U:0"])
vertLocList.append(df["Points:2"]*100)
print(vertLocList[0])
horLocList.append(df["U:0"]*500+df["Points:0"]*100)
print(horLocList[0])



figure(figsize=(20, 4), dpi=300)
plt.plot(x,z, marker="o", linestyle = "None", markersize=2)
plt.plot(horLocList[0],vertLocList[0], marker="o", linestyle="None", markersize=2)
plt.ylim(0,16)
plt.xlim(0,760)
plt.savefig("testPlot.png")

