from pysb_Boolean import *
from pysb.importers.boolean import model_from_boolean

mode = 'ROA'  # 'GSP', 'GA', 'ROA'
model = model_from_boolean('mapk_soce_oscillate.txt', mode=mode)
n_runs = 1000

step_labels = [
    "equilibration",
    "add BRAFi"
]
delta_ts = [50, 500]
conditions = [
    None,
    [("BRAFi", True)]
]

sim_steps = [SimStep(label, delta_t, condition) for label, delta_t, condition
             in zip(step_labels, delta_ts, conditions)]

tspans, outputs = sim_protocol(model, sim_steps, n_runs=n_runs, t_start=-sim_steps[0].delta_t, mode=mode,
                               verbose=True)

obs_names = [
    ['Ca_cyt_1_True_obs', 'Ca_cyt_2_True_obs', 'Ca_cyt_3_True_obs'],
    'BRAF_True_obs',
    'Ca_channel_True_obs',
    'Ca_ER_True_obs',
    'Ca_ext_True_obs',
    'Ca_pump_ER_True_obs',
    'ERK_True_obs',
    'Gene_exp_True_obs',
    'MEK_True_obs'
]
obs_labels = ['Ca_cyt', 'BRAF', 'Ca_channel', 'Ca_ER', 'Ca_ext', 'Ca_pump_ER', 'ERK', 'Gene_exp', 'MEK']
obs_colors = ['blue', 'green', 'black', 'red', 'purple', 'brown', 'yellow', 'orange', 'cyan']
# obs_markers = ['s', 'o', '^', '*', 'd', 'H', 'v', '<', '>']

observables = [ObsToPlot(name, label, color) for name, label, color in zip(obs_names, obs_labels, obs_colors)]

plot_results(tspans, outputs, observables, mode, multi_plots=False, show_plots=True)
