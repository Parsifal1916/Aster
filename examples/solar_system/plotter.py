import pandas as pd
import matplotlib.pyplot as plt

# Carica i file CSV
x_data = pd.read_csv("Graph2.csv")
y_data = pd.read_csv("Graph3.csv")

# Controlla che i corpi siano gli stessi
x_bodies = [col for col in x_data.columns if col != 'time']
y_bodies = [col for col in y_data.columns if col != 'time']

assert x_bodies == y_bodies, "I corpi nei due file non corrispondono."
bodies = x_bodies

# Tronca al numero minimo di righe valide comuni
min_len = min(len(x_data), len(y_data))
x_data = x_data.iloc[:min_len]
y_data = y_data.iloc[:min_len]

# Colori distinti per ogni corpo (più se servono)
color_palette = ['red', 'green', 'blue', 'orange', 'purple', 'cyan', 'magenta']

# Plot
plt.figure(figsize=(8, 6))

for i, body in enumerate(bodies):
    color = color_palette[i % len(color_palette)]
    plt.plot(x_data[body], y_data[body], label=body, color=color)

plt.xlabel('x')
plt.ylabel('y')
plt.title('Traiettorie dei corpi')
plt.legend()
plt.grid(True)
plt.axis('equal')  # Mantiene le proporzioni
plt.tight_layout()
plt.show()
