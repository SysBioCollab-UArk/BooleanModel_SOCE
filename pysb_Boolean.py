from pysb.importers.boolean import model_from_boolean
from pysb.simulator import BngSimulator
import numpy as np
import matplotlib.pyplot as plt
from pysb.bng import *
from pysb.generator.bng import BngPrinter
from pysb.bng import BngFileInterface

mode = 'GSP'  # 'ROA'
model = model_from_boolean('mapk_soce.txt', mode=mode)

n_runs = 1000

species_to_plot = ['Ca_cyt', 'BRAF', 'Ca_channel', 'Ca_ER', 'Ca_ext', 'Ca_pump_ER', 'ERK', 'Gene_exp', 'MEK']
observables = [
    ['Ca_cyt_1_True_obs', 'Ca_cyt_2_True_obs', 'Ca_cyt_3_True_obs'],
    ['BRAF_True_obs'],
    ['Ca_channel_True_obs'],
    ['Ca_ER_True_obs'],
    ['Ca_ext_True_obs'],
    ['Ca_pump_ER_True_obs'],
    ['ERK_True_obs'],
    ['Gene_exp_True_obs'],
    ['MEK_True_obs']
]
markers = ['s', 'o', '^', '*', 'd', 'H', 'v', '<', '>']
colors = ['blue', 'green', 'black', 'red', 'purple', 'brown', 'yellow', 'orange', 'cyan']

# tspan1 = np.linspace(0, 10, 101)
# sim = BngSimulator(model, verbose=True)
# out = sim.run(tspan=tspan1, n_runs=2)
# print(np.array(out.observables).shape)
# quit()

# time intervals for each step
tend_1 = 10  # step1: equilibration
tend_2 = 50  # step2: add BRAFi
tend_3 = 10  # step3: remove external calcium and add pump inhibitor
tend_4 = 50  # step4: add external calcium
tend_5 = 50  # step5: remove pump inhibitor

with BngFileInterface(model, verbose=True) as bng:
    bng.action('generate_network')
    for n in range(n_runs):
        # print(n)
        # step1: equilibration
        tspan1 = np.linspace(0, tend_1, tend_1+1)
        tspan_all = tspan1
        bng.action('simulate', method='ssa', t_start=tspan_all[0], t_end=tspan_all[-1],
                   n_steps=len(tspan1)-1, suffix=str(n))

        # step2: add BRAFi
        tspan2 = np.linspace(0, tend_2, tend_2+1)
        bng.set_concentration("BRAFi(state~False)", 0)
        bng.set_concentration("BRAFi(state~True)", 1)
        kwarg = {'continue': 1}
        bng.action('simulate', method='ssa', t_start=tspan_all[-1], t_end=tspan_all[-1]+tspan2[-1],
                   n_steps=len(tspan2)-1, suffix=str(n), **kwarg)
        tspan_all = np.append(tspan_all, tspan2[1:] + tspan1[-1])

        # step3: remove external calcium and add pump inhibitor
        tspan3 = np.linspace(0, tend_3, tend_3+1)
        bng.set_concentration("Ca_ext(state~False)", 1)
        bng.set_concentration("Ca_ext(state~True)", 0)
        bng.set_concentration("pumpi(state~False)", 0)
        bng.set_concentration("pumpi(state~True)", 1)
        bng.action('simulate', method='ssa', t_start=tspan_all[-1], t_end=tspan_all[-1]+tspan3[-1],
                   n_steps=len(tspan3)-1, suffix=str(n), **kwarg)
        tspan_all = np.append(tspan_all, tspan3[1:] + tspan1[-1] + tspan2[-1])

        # step4: add external calcium
        tspan4 = np.linspace(0, tend_4, tend_4+1)
        bng.set_concentration("Ca_ext(state~False)", 0)
        bng.set_concentration("Ca_ext(state~True)", 1)
        bng.action('simulate', method='ssa', t_start=tspan_all[-1], t_end=tspan_all[-1]+tspan4[-1],
                   n_steps=len(tspan4)-1, suffix=str(n), **kwarg)
        tspan_all = np.append(tspan_all, tspan4[1:] + tspan1[-1] + tspan2[-1] + tspan3[-1])

        # step5: remove pump inhibitor
        tspan5 = np.linspace(0, tend_5, tend_5+1)
        bng.set_concentration("pumpi(state~False)", 1)
        bng.set_concentration("pumpi(state~True)", 0)
        bng.action('simulate', method='ssa', t_start=tspan_all[-1], t_end=tspan_all[-1]+tspan5[-1],
                   n_steps=len(tspan5)-1, suffix=str(n), **kwarg)
        tspan_all = np.append(tspan_all, tspan5[1:] + tspan1[-1] + tspan2[-1] + tspan3[-1] + tspan4[-1])

        # reset concentration
        bng.action('resetConcentrations')

    # run the simulation and plot the results
    bng.execute()
    result = BngFileInterface.read_simulation_results_multi(['%s_%d' % (bng.base_filename, n) for n in range(n_runs)])
    result = np.array(result)
    plt.figure(figsize=(12.8, 4.8))
    for sp, obs, marker, color in zip(species_to_plot, observables, markers, colors):
        print(sp)
        avg_traj = np.mean(result[:][obs[0]], axis=0)
        # n_obs = len(obs)
        if len(obs) > 1:
            for i in range(1, len(obs)):
                avg_traj += np.mean(result[:][obs[i]], axis=0)
            avg_traj /= len(obs)
        # for n in range(n_runs):
        #     plt.plot(tspan_all, result[n][obs[0]], lw=1, color=obs[1])
        plt.plot(tspan_all[len(tspan1):]-tspan1[-1], avg_traj[len(tspan1):], marker=marker, color=color, label=sp)

plt.legend(loc=0, fontsize=15)
plt.xlim(xmin=0)
plt.ylim((-0.1, 1.1))
plt.xlabel('time', fontsize=16)
plt.ylabel('probability ON', fontsize=16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.tight_layout()

plt.show()
