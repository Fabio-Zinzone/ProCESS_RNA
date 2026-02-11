import tellurium as te
import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np

output_folder = "DEBUG"
sbml_path = os.path.join(output_folder, "model.xml")
results_path = os.path.join(output_folder, "final.csv")
data_path = os.path.join(output_folder, "measurement.csv")

if not os.path.exists(sbml_path):
    raise FileNotFoundError("Missing model.xml!")
r = te.loadSBMLModel(sbml_path)

df_res = pd.read_csv(results_path)
df_data = pd.read_csv(data_path)

print(">>> Parameter insertion...")
for _, row in df_res.iterrows():
    pid = row["Parameter"]
    val = row["Estimated Value"]

    if not np.isnan(val) and pid in r.model:
        r.model[pid] = val

print(">>> Best fit simulation...")
sim = r.simulate(0, 100, 1000)
df_sim = pd.DataFrame(sim, columns=sim.colnames)
df_sim.columns = [c.replace("[", "").replace("]", "") for c in df_sim.columns]

print(">>> Graph generation...")
genes = [f"G{i}" for i in range(20)]
cols = 4
rows = 5
fig, axes = plt.subplots(rows, cols, figsize=(20, 15))
fig.suptitle("Data vs Fit", fontsize=16)
axes = axes.flatten()

for i, gene in enumerate(genes):
    ax = axes[i]
    
    if gene in df_data.columns:
        ax.scatter(df_data["time"][::15], df_data[gene][::15], 
                   color="red", alpha=0.5, s=20, label="Data" if i==0 else "")
    
    if gene in df_sim.columns:
        ax.plot(df_sim["time"], df_sim[gene], 
                color="blue", linewidth=2, label="Fit" if i==0 else "")
    
    k_param = f"K{i}"
    err_str = ""
    try:
        row = df_res[df_res["Parameter"] == k_param].iloc[0]
        err = ((row['Estimated Value'] - row['True Value']) / row['True Value']) * 100
        err_str = f"Err K: {err:.1f}%"
    except: pass

    ax.set_title(f"{gene} ({err_str})", fontsize=10, fontweight="bold")
    ax.grid(True, alpha=0.3)

    if i == 0: ax.legend()

plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig(os.path.join(output_folder, "fit.png"))
print(f"\n>>> GRAPH SAVED: {os.path.join(output_folder, 'fit.png')}")
