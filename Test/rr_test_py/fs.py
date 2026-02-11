import matplotlib.pyplot as plt
import numpy as np
import io
import pandas as pd

data_raw_25 = """
2    6.364158    6.691285    6.170957    9.028159    12.443324    13.214243    0.356019    1.887862    0.091792    0.407260
51   6.364150    6.691266    6.170947    9.028128    12.443306    13.214201    0.356021    1.887865    0.091792    0.407260 
1000 6.364150    6.691267    6.170949    9.028135    12.443308    13.214208    0.356021    1.887865    0.091792    0.407260
10000 6.364155   6.691277    6.170953    9.028147    12.443316    13.214227    0.356020    1.887863    0.091792    0.407260
"""

data_raw_5000 = """
2    10.099020    40.396078    10.098468    50.492341    20.097488    100.487439    0.000245    0.000490    0.100000    0.800000
51   10.099020    40.396078    10.098468    50.492341    20.097488    100.487439    0.000245    0.000490    0.100000    0.800000 
1000 10.099020    40.396078    10.098468    50.492341    20.097488    100.487439    0.000245    0.000490    0.100000    0.800000
10000 10.099020   40.396078    10.098468    50.492341    20.097488    100.487439    0.000245    0.000490    0.100000    0.800000
"""

species_names = ['m_E1', 'p_E1', 'm_COL1a', 'p_COL1a', 'm_FT4', 'p_FT4', 'm_FT2a', 'p_FT2a', 'm_AP1a', 'p_AP1a']
step_labels = {
    2: '2 Steps',
    51: 'Default (51)',
    1000: '1000 Steps',
    10000: '10000 Steps'
}

df = pd.read_csv(io.StringIO(data_raw_25), sep='\s+', engine='python', header=None)
df.columns = ['N_Steps'] + species_names
df = df.set_index('N_Steps')

bar_width = 0.25
fig, ax = plt.subplots(figsize=(14, 7))

# Iterazione sull'indice numerico del DataFrame
for i, idx in enumerate(df.index):
    step_label = step_labels[idx]
    x_positions = np.arange(len(species_names)) + (i * bar_width)
    values = df.loc[idx] 
    ax.bar(x_positions, values, bar_width, label=step_label)

ax.set_xticks(np.arange(len(species_names)) + bar_width * 1.5)
ax.set_xticklabels(species_names, rotation=45, ha='right')
ax.set_title('Final State Comparison (t=25)')
ax.set_ylabel('Amounts')
ax.legend()
ax.grid(axis='y', linestyle='--', alpha=0.6)

plt.tight_layout()
plt.savefig('final_state.png')
