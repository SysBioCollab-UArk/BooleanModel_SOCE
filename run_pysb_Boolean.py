import numpy as np
from pysb.importers.boolean import model_from_boolean
from pysb_Boolean import *
from util import get_sim_steps
import os

path = 'VERSIONS'
files = sorted([file for file in os.listdir(path) if os.path.isfile(os.path.join(path, file))])

for file in ['mapk_soce_v6.txt']:  # files:

    print('model: %s' % file)

    mode = 'GSP'  # 'GSP', 'GA', 'ROA'
    model = model_from_boolean(os.path.join(path, file), mode=mode)
    n_runs = 1000

    # conditions = ["untreated", "30_min_BRAFi", "8_days_BRAFi", "8_days_BRAFi_plus_MEKi"]
    # xlims = [(40, 70), (40, 70), (290, 320), (300, 330)]

    conditions = ["8_days_BRAFi"]
    xlims = [(0, 30)]

    for condition, xlim in zip(conditions, xlims):

        print('=== %s ===' % condition)

        sim_steps = get_sim_steps(condition)

        # Loop over a range of values for when ER pump inhibitor is added
        for dt in np.arange(250, 350+1, 25):
            print('** delta_t =', dt)
            sim_steps[1].delta_t = dt
            tspans, outputs = sim_protocol(model, sim_steps, n_runs=n_runs,
                                           t_start=-sim_steps[0].delta_t-sim_steps[1].delta_t+10,
                                           mode=mode, verbose=False)

            obs_names = [
                ['Ca_cyt_1_True_obs', 'Ca_cyt_2_True_obs', 'Ca_cyt_3_True_obs'] # ,
                # 'BRAF_True_obs',
                # 'Ca_channel_True_obs',
                # 'Ca_ER_True_obs',
                # 'Ca_ext_True_obs',
                # 'Ca_pump_ER_True_obs',
                # 'ERK_True_obs',
                # 'Gene_exp_True_obs',
                # 'MEK_True_obs'
            ]
            obs_labels = ['Ca_cyt']  # , 'BRAF', 'Ca_channel', 'Ca_ER', 'Ca_ext', 'Ca_pump_ER', 'ERK', 'Gene_exp', 'MEK']
            obs_colors = ['blue']  # , 'green', 'black', 'red', 'purple', 'brown', 'yellow', 'orange', 'cyan']
            # obs_markers = ['s', 'o', '^', '*', 'd', 'H', 'v', '<', '>']

            observables = [ObsToPlot(name, label, color) for name, label, color in
                           zip(obs_names, obs_labels, obs_colors)]

            outfile_name = '%s_%s.pdf' % (file[:-4], condition)  # remove the .txt extension from the filename
            plot_results(tspans, outputs, observables, multi_plots=False, save_plots=outfile_name, show_plots=False,
                         figname='%s_%s' % (file[:-4], condition), xlim=xlim, ylim=(0.15, 1.05), lw=1)
