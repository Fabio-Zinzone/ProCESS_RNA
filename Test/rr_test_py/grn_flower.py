import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import tellurium as te
import roadrunner
import time

start = time.time()

g1, g2, g3, g4, g5 = 'E1', 'COL1a', 'FT4', 'FT2a', 'AP1a'
cg3 = "[FT4]"
m, p = "p_", "m_" "[m_E1]"
g = [g1, g2, g3, g4, g5]
c = ['g', 'r', 'b', 'c', 'm']
lw = 1
sn = 0.1
ni, nm, nf, ns = 0, 25, 50, 1000
intg = 'cvode' # cvode gillespie rk4 rk45

mpl.rcParams.update({
    'font.size': 24,
    'axes.titlesize': 34,
    'axes.labelsize': 26,
    'legend.fontsize': 22,
    'xtick.labelsize': 22,
    'ytick.labelsize': 22
})

file_antimony = "grn_flower.txt"

with open(file_antimony, 'r') as f:
        antimony_str = f.read()

sbml_str = te.antimonyToSBML(antimony_str)
file_xml = "grn_flower.xml"
with open(file_xml, "w") as f:
    f.write(sbml_str)

col_list = ['time']
for name in g:
    col_list.append(f'm_{name}')
    col_list.append(f'p_{name}')

rr = roadrunner.RoadRunner("grn_flower.xml")
rr.setIntegrator(intg)
rr.timeCourseSelections = col_list
rr.structured_results = True
startx = time.time()
ts = rr.simulate(ni,nm,ns)
stopx = time.time()
nome_file = "time_series.txt"
np.savetxt(nome_file, ts, fmt='%.6f', delimiter='\t')
df = pd.DataFrame(ts, columns=rr.timeCourseSelections).set_index('time')
print(df)

sbml_str_agg = rr.getCurrentSBML()
with open("grn_flower_agg.xml", "w") as f:
    f.write(sbml_str_agg)

antimony_str_agg = te.sbmlToAntimony(sbml_str_agg)
with open("grn_flower_agg.txt", "w") as f:
    f.write(antimony_str_agg)

n = np.random.normal(0, sn, size=ts.shape)
n[:, 0] = 0
ts += n # *= (1.0 + n)

rrx = roadrunner.RoadRunner("grn_flower_agg.xml")
rrx.setIntegrator(intg)
rrx.timeCourseSelections = col_list
rrx.structured_results = False
tsx = rrx.simulate(nm,nf,ns)
nome_file = "time_series_agg.txt"
np.savetxt(nome_file, tsx, fmt='%.6f', delimiter='\t')
dfx = pd.DataFrame(tsx, columns=rrx.timeCourseSelections).set_index('time')
print(dfx)

nx = np.random.normal(0, sn, size=tsx.shape)
nx[:, 0] = 0
tsx += n # *= (1.0 + nx)

plt.figure(figsize=(36,18)) # mRNA 2 4 6 8 10 - 1 3 5 7 9
for i in range(5):
    ci = 2 * i + 1
    plt.plot(ts[:, 0], ts[:, ci], label=g[i], color=c[i], linewidth=lw)
plt.title('First Run Time Series: mRNA')
plt.xlabel('t')
plt.ylabel('mRNA')
plt.legend()
plt.grid()

plt.figure(figsize=(36,18)) # protein 3 5 7 9 11 - 2 4 6 8 10
for i in range(5):
    ci = 2 * i + 2
    plt.plot(ts[:, 0], ts[:, ci], label=g[i], color=c[i], linewidth=lw)
plt.title('First Run Time Series: Protein') 
plt.xlabel('t')
plt.ylabel('Protein')
plt.legend()
plt.grid()

plt.figure(figsize=(36,18)) # mRNA 2 4 6 8 10 - 1 3 5 7 9
for i in range(5):
    ci = 2 * i + 1
    plt.plot(tsx[:, 0], tsx[:, ci], label=g[i], color=c[i], linewidth=lw)
plt.title('Second Run Time Series: mRNA')
plt.xlabel('t')
plt.ylabel('mRNA')
plt.legend()
plt.grid()

plt.figure(figsize=(36,18)) # protein 3 5 7 9 11 - 2 4 6 8 10
for i in range(5):
    ci = 2 * i + 2
    plt.plot(tsx[:, 0], tsx[:, ci], label=g[i], color=c[i], linewidth=lw)
plt.title('Second Run Time Series: Protein') 
plt.xlabel('t')
plt.ylabel('Protein')
plt.legend()
plt.grid()

with PdfPages('grn_flower_'+intg+'.pdf') as pdf:
    for fig_num in plt.get_fignums():
        pdf.savefig(fig_num)

stop = time.time()
delta_t = stop - start
delta_tx = stopx - startx
print("Step number: ", ns)
print("Execution time: ", delta_t)
print("Execution time simulation: ", delta_tx)

#print(rr.getFloatingSpeciesIds())
#print(rr.getGlobalParameterIds())

#print(rr["m_"+g3])
#print(rr["[m_"+g3+"]"])

#plt.show()
