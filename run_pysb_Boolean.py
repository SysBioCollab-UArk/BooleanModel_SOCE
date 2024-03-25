from pysb.importers.boolean import model_from_boolean
from pysb_Boolean import *
from util import get_sim_steps

versions = [1, 2, 3, 4]

for version in versions:

    print('version: %d' % version)

    mode = 'GSP'  # 'GSP', 'GA', 'ROA'
    model = model_from_boolean('VERSIONS/mapk_soce_V%d.txt' % version, mode=mode)
    n_runs = 100

    conditions = ["untreated", "30_min_BRAFi", "8_days_BRAFi", "8_days_BRAFi_plus_MEKi"]
    xlims = [(40, 70), (40, 70), (490, 520), (500, 530)]

    for condition, xlim in zip(conditions, xlims):

        print('=== %s ===' % condition)

        sim_steps = get_sim_steps(condition)
        tspans, outputs = sim_protocol(model, sim_steps, n_runs=n_runs, t_start=-sim_steps[0].delta_t, mode=mode,
                                       verbose=False)

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

        observables = [ObsToPlot(name, label, color) for name, label, color in zip(obs_names, obs_labels, obs_colors)]

        outfile_name = 'V%d_%s.pdf' % (version, condition)
        plot_results(tspans, outputs, observables, mode, multi_plots=False, save_plots=outfile_name, show_plots=False,
                     xlim=xlim, ylim=(0.15, 1.05))