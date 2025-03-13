import pandas as pd
import matplotlib.pyplot as plt

x_df = pd.read_csv("Graph1.csv")
y_df = pd.read_csv("Graph2.csv")

body_labels = x_df.columns[1:] 

plt.figure(figsize=(8, 8))

for body in body_labels:
    plt.plot(x_df[body], y_df[body])

plt.xlabel("X")
plt.ylabel("Y")
plt.title("3-body PN1 ")
plt.grid()
plt.savefig("3body_PN1.png", dpi=300, bbox_inches="tight")
plt.show()