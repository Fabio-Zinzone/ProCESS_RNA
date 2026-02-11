import os
import sys
import time
import tellurium as te
import pandas as pd
import numpy as np
import petab
import pypesto.petab
import pypesto.optimize
import amici
import logging
import warnings
import traceback

start = time.time()

np.random.seed(0)

logging.getLogger("pypesto").setLevel(logging.WARNING)
logging.getLogger("amici").setLevel(logging.ERROR)
logging.getLogger("petab").setLevel(logging.ERROR)
warnings.filterwarnings("ignore")

print("\n" + "="*60)
print("   AMICI TEST   ")
print("="*60 + "\n")

output_folder = "DEBUG"
os.makedirs(output_folder, exist_ok=True)
model_file = "anti_20.txt"                                                          # antimony test string
n_genes = 20

if not os.path.exists(model_file):                                                  # model
    print(f"ERROR: File '{model_file}' not found.")
    sys.exit(1)

with open(model_file, "r") as f:
    antimony_str = f.read()

rr = te.loada(antimony_str)                                                         # antimony conversion with tellurium

true_params = {}
all_params = rr.model.getGlobalParameterIds()

for pid in all_params:
    true_params[pid] = rr.model[pid]                                                # store true values

for pid in all_params:
    if pid == "n" or pid == "lam":
        rr.model[pid] = true_params[pid]                                            # fixed values
    else:
        rr.model[pid] = 1.0                                                         # initial condition for sbml

sbml_path = os.path.join(output_folder, "model.xml")                             # model input
with open(sbml_path, "w") as f:
    f.write(rr.getSBML())

for pid, val in true_params.items():
    if pid in rr.model: rr.model[pid] = val                                  # restore params for simulation

print(">>> Synthetic data generation...")
data_true = rr.simulate(0, 100, 1000)                                               # roadrunner simulation

df_noisy = pd.DataFrame(data_true, columns=data_true.colnames)                      # noise dataframe
df_noisy.columns = [c.replace("[", "").replace("]", "") for c in df_noisy.columns]
noise_cols = [c for c in df_noisy.columns if c != "time"]

df_noisy[noise_cols] += np.random.normal(0, 0.005, df_noisy[noise_cols].shape)      # gaussian noise
df_noisy[df_noisy < 0] = 0.001
df_noisy.to_csv(os.path.join(output_folder, "measurement.csv"), index=False)        # save measurement file

pd.DataFrame(list(true_params.items()), columns=["Parameter", "True Value"]).to_csv(
    os.path.join(output_folder, "real_par.csv"), index=False
)

# PETAB CONFIGURATION

df_measurement = df_noisy.melt(id_vars="time", var_name="observableId", value_name="measurement")
df_measurement["simulationConditionId"] = "cond1"
df_measurement["observableId"] = "obs_" + df_measurement["observableId"]

# observables input

observables = [{
    "observableId": f"obs_G{i}",
    "observableFormula": f"G{i}",
    "noiseFormula": "0.005"
} for i in range(n_genes)]
df_observable = pd.DataFrame(observables).set_index("observableId")
df_observable.to_csv(os.path.join(output_folder, "observable.csv"), index=True)      # save observables

# parameter input

parameters_petab = []
for pid in all_params:

    if pid == "n" or pid == "lam":
        estimate = 0
        nominal_value = true_params[pid]
        scale = "lin"                                                                  # linear scale: n, lam
        lb, ub = 0.1, 10.0
        
    else:
        estimate = 1
        nominal_value = 1.0
        
        if pid.startswith("h"):
            scale = "log10"                                                             # log scale: h
            lb, ub = 0.1, 100.0
        else: 
            scale = "lin"                                                               # linear scale: K
            lb, ub = 0.1, 15.0

    parameters_petab.append({
        "parameterId": pid,
        "parameterScale": scale,
        "lowerBound": lb,
        "upperBound": ub,
        "estimate": estimate, 
        "nominalValue": nominal_value 
    })

df_parameter = pd.DataFrame(parameters_petab).set_index("parameterId")
df_parameter.to_csv(os.path.join(output_folder, "parameter.csv"), index=True)          # save parameter file

# condition input

df_condition = pd.DataFrame({"conditionId": ["cond1"]}).set_index("conditionId")
df_condition.to_csv(os.path.join(output_folder, "condition.csv"), index=True)          # save condition file

petab_problem = petab.Problem(                                                          # PETAB problem
    model=petab.models.sbml_model.SbmlModel.from_file(sbml_path),
    condition_df=df_condition,
    measurement_df=df_measurement,
    parameter_df=df_parameter,
    observable_df=df_observable
)

importer = pypesto.petab.PetabImporter(petab_problem)                                   # import PyPESTO                         
model = importer.create_model(force_compile=True) 
solver = model.getSolver()
solver.setMaxSteps(20000)
solver.setSensitivityMethod(amici.SensitivityMethod_forward)                            # forward sensitivity
#solver.setSensitivityMethod(amici.SensitivityMethod_adjoint)
solver.setMaxTime(30.0)
solver.setAbsoluteTolerance(1e-4)
solver.setRelativeTolerance(1e-4)

problem = importer.create_problem(model=model, solver=solver)

print(">>> SELF-CHECK...")                                                              # self-check
try:
    x_start = []
    for pid in petab_problem.parameter_df.index:
        if petab_problem.parameter_df.loc[pid, "estimate"] == 1:
            val = petab_problem.parameter_df.loc[pid, "nominalValue"]
            scale = petab_problem.parameter_df.loc[pid, "parameterScale"]
            x_start.append(np.log10(val) if scale == "log10" else val)
    
    x_start = np.array(x_start)
    print(f"   Parameters: {len(x_start)}")
    
    llh = problem.objective(x_start)
    print(f"   SELF-CHECK: Log-Likelihood = {llh}")
    
    if np.isnan(llh) or np.isinf(llh):
        print("   FATAL ERROR: Simulation failure (NaN/Inf).")
        sys.exit(1)
        
except Exception as e:
    print(f"   ERROR DURING SELF-CHECK: {e}")
    traceback.print_exc()
    sys.exit(1)

print(">>> SELF-CHECK FINISHED. Begin optimization...")

optimizer = pypesto.optimize.ScipyOptimizer(                                            # PyPESTO optimizer
    method="l-bfgs-b",
    options={
        'maxiter': 1000,                                                                # max iter
        'ftol': 1e-8
    }
)

result = pypesto.optimize.minimize(                                                     # PyPESTO optimization
    problem=problem,
    optimizer=optimizer,
    n_starts=20,                                                                        # independent starts
    filename=None,
    progress_bar=True
)

best_result = result.optimize_result.list[0]                                            # get best result
final_data = []

for pid in true_params.keys():
    val_true = true_params[pid]
    is_estimated = df_parameter.loc[pid, "estimate"] == 1
    val_est = np.nan
    
    if is_estimated:
        try:
            if pid in problem.x_names:
                idx = problem.x_names.index(pid)
                raw_val = best_result['x'][idx]
                if df_parameter.loc[pid, "parameterScale"] == "log10":
                    val_est = 10**raw_val
                else:
                    val_est = raw_val
        except: pass
    else:
        val_est = df_parameter.loc[pid, "nominalValue"]

    final_data.append({
        "Parameter": pid,
        "True value": val_true,
        "Estimated value": val_est,
        "Absolute error": abs(val_est - val_true) if not np.isnan(val_est) else np.nan,
        "Type": "Estimated" if is_estimated else "Fixed"
    })

df_final = pd.DataFrame(final_data)
df_final.to_csv(os.path.join(output_folder, "final.csv"), index=False)
print(f"\n   ANALYSIS COMPLETED. Data saved in '{output_folder}/final.csv'.")           # "final" saving
print()

stop = time.time()
delta_t = stop - start
print(f"Execution time: {delta_t} s = {delta_t/60} m = {delta_t/3600} h")
print()
