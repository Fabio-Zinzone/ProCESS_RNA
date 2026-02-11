import matplotlib.pyplot as plt
import numpy as np

n_steps = [2, 10, 100, 1000, 10000, 100000, 1000000]
t_sim   = [0.000592, 0.000627, 0.001130, 0.004429, 0.039267, 0.367768, 3.521741]
t_tot   = [1.162974, 1.193904, 1.289451, 1.366338, 2.430910, 9.455502, 49.471379]

plt.figure(figsize=(10, 6))
plt.semilogx(n_steps, t_tot, 'r-o', label='Tempo Totale (Programma)', linewidth=2)
for x, y in zip(n_steps, t_tot):
    plt.text(x, y, f"{y:.2f}s\n[x{y/t_tot[0]:.1f}]", fontsize=10, ha='right', va='bottom', color='darkred')
plt.xlabel('Numero Step')
plt.ylabel('Tempo')
plt.legend()
plt.grid(True, which="both", ls="-", alpha=0.5)

plt.figure(figsize=(10, 6))
plt.semilogx(n_steps, t_sim, 'b-s', label='Tempo Simulazione (RoadRunner)', linewidth=2)
for x, y in zip(n_steps, t_sim):
    plt.text(x, y, f"{y:.6f}s\n[x{y/t_sim[0]:.1f}]", fontsize=10, ha='right', va='bottom', color='darkblue')
plt.xlabel('Numero Step')
plt.ylabel('Tempo')
plt.legend()
plt.grid(True, which="both", ls="-", alpha=0.5)

plt.tight_layout()
plt.show()
