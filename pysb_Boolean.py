from pysb.importers.boolean import model_from_boolean
from pysb.simulator import BngSimulator
import numpy as np
import matplotlib.pyplot as plt


class SimStep(object):
    def __init__(self, label, delta_t, condition=None):
        self.label = label
        self.delta_t = delta_t
        self.condition = condition


def sim_protocol(model, sim_steps, n_runs=1, t_start=0, mode='GSP', verbose=False):
    sim = BngSimulator(model, verbose=verbose)
    tspans = []
    outputs = []
    additional_args = {}
    end_string = ""
    sp_names = None
    if mode == 'GSP' or mode == 'GA':
        end_string = ")"
    elif mode == 'ROA':
        end_string = ", reset='N')"
        # The output_step_interval is the number of nodes with update rules, multiplied by 2. That includes the model
        # monomers (not necessarily all of them) PLUS the RESET monomer that's created automatically.
        additional_args['output_step_interval'] = (int(model.parameters['N_NODES'].value) + 1) * 2
    else:
        print("Mode '%s' not recognized. Quitting.")
        quit(1)
    for step, ss in enumerate(sim_steps):
        print("step %d: %s" % (step, ss.label))
        # first step
        if step == 0:
            initials = None
            if ss.condition is not None:
                initials = {}
                for c in ss.condition:
                    if c is not None:
                        initials["%s(state='%s'%s" % (c[0], str(c[1]), end_string)] = 1
                        initials["%s(state='%s'%s" % (c[0], str(not c[1]), end_string)] = 0
        # subsequent steps
        else:
            t_start = tspans[-1][-1]
            initials = [sp[-1] for sp in outputs[-1].species]
            n_runs = 1
            if ss.condition is not None:
                for c in ss.condition:
                    idx1 = sp_names.index("%s(state='%s'%s" % (c[0], str(c[1]), end_string))
                    idx2 = sp_names.index("%s(state='%s'%s" % (c[0], str(not c[1]), end_string))
                    for init in initials:
                        init[idx1] = 1
                        init[idx2] = 0
        # run simulations
        if mode == 'GSP' or mode == 'GA':
            tspans.append(np.linspace(t_start, t_start + ss.delta_t, int(ss.delta_t * 10) + 1))
            time = tspans[-1]
        elif mode == 'ROA':
            additional_args['max_sim_steps'] = additional_args['output_step_interval'] * ss.delta_t
            tspans.append(t_start + np.arange(ss.delta_t + 1))
            time = [0, additional_args['max_sim_steps'] * 1000]
        outputs.append(sim.run(tspan=time, initials=initials, n_runs=n_runs, cleanup=True, **additional_args))
        # store species names as strings if this is the first simulation
        if sp_names is None:
            sp_names = [str(m) for m in model.species]

    return tspans, outputs


def get_sim_steps(which):

    _supported = ["untreated", "30_min_BRAFi", "8_days_BRAFi", "8_days_BRAFi_plus_MEKi"]

    if which not in _supported:
        print("Error in 'pysb_Boolean::get_sim_steps': requested protocol ('%s') not found. Please try again." %
              str(which))
        quit()

    step_labels = None
    delta_ts = None
    conditions = None

    # Untreated
    if which == "untreated":

        step_labels = [
            "equilibration",
            "remove external calcium",
            "add pump inhibitor",
            "add external calcium"
        ]
        delta_ts = [50, 50, 10, 100]
        conditions = [
            None,
            [("Ca_ext", False)],
            [("pumpi", True)],
            [("Ca_ext", True)]
        ]

    # 30 min BRAFi
    elif which == "30_min_BRAFi":

        step_labels = [
            "equilibration",
            "remove external calcium and add BRAFi",
            "add pump inhibitor",
            "add external calcium"
        ]
        delta_ts = [50, 50, 10, 100]
        conditions = [
            None,
            [("Ca_ext", False), ("BRAFi", True)],
            [("pumpi", True)],
            [("Ca_ext", True)]
        ]

    # 8 days BRAFi
    elif which == "8_days_BRAFi":

        step_labels = [
            "equilibration",
            "remove external calcium and add BRAFi",
            "add pump inhibitor",
            "add external calcium"
        ]
        delta_ts = [50, 500, 10, 100]
        conditions = [
            None,
            [("Ca_ext", False), ("BRAFi", True)],
            [("pumpi", True)],
            [("Ca_ext", True)]
        ]

    # 8 days BRAFi + 15 min MEKi
    elif which == "8_days_BRAFi_plus_MEKi":

        step_labels = [
            "equilibration",
            "remove external calcium and add BRAFi",
            "add MEKi",
            "add pump inhibitor",
            "add external calcium"
        ]
        delta_ts = [50, 500, 10, 10, 100]
        conditions = [
            None,
            [("Ca_ext", False), ("BRAFi", True)],
            [("MEKi", True)],
            [("pumpi", True)],
            [("Ca_ext", True)]
        ]

    return [SimStep(label, dt, cond) for label, dt, cond in zip(step_labels, delta_ts, conditions)]


class ObsToPlot(object):
    def __init__(self, name, label, color=None, marker=None):
        if isinstance(name, str):
            name = [name]
        self.name = name
        self.label = label
        self.color = color
        self.marker = marker


def plot_results(tspans, outputs, observables, mode, multi_plots=False, show_plots=False):
    if not multi_plots:
        plt.figure(figsize=(12.8, 4.8))
        plt.ylim(bottom=-0.05, top=1.05)
    lines = []
    labels = []
    for obs in observables:
        if multi_plots:
            plt.figure(figsize=(12.8, 4.8))
            plt.ylim(bottom=-0.05, top=1.05)
            lines = []
            labels = []
        labels.append(obs.label)
        for step, out, tspan in zip(range(len(outputs)), outputs, tspans):
            out_observables = np.array(out.observables)
            y = np.mean(out_observables[obs.name[0]], axis=0)
            if len(obs.name) > 1:
                for i in range(1, len(obs.name)):
                    y += np.mean(out_observables[obs.name[i]], axis=0)
                y /= len(obs.name)
            line, = plt.plot(tspan, y, lw=3, color=obs.color)
            if step == 0:
                lines.append(line)
            # plot vertical lines for each step
            bottom, top = plt.ylim()
            plt.plot([tspan[-1], tspan[-1]], [bottom, top], 'k--')
    plt.xlabel('round' if mode == 'ROA' else 'time', fontsize=16)
    plt.ylabel('value', fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.legend(loc=0, handles=lines, labels=labels, bbox_to_anchor=(1, 1), fontsize=15)
    plt.tight_layout()
    if show_plots:
        plt.show()


if __name__ == '__main__':

    mode = 'GSP'  # 'GSP', 'GA', 'ROA'
    model = model_from_boolean('VERSIONS/mapk_soce_V1.txt', mode=mode)
    n_runs = 100

    step_labels = [
        "equilibration",
        "add BRAFi",
        "remove external calcium and add pump inhibitor",
        "add external calcium"  # ,
        # "remove pump inhibitor"
    ]
    delta_ts = [40, 100, 10, 100]
    conditions = [
        None,
        [("BRAFi", True)],
        [("Ca_ext", False), ("pumpi", True)],
        [("Ca_ext", True)],  # ,
        # [("pumpi", False)],
        [["MEKi", True]]
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
