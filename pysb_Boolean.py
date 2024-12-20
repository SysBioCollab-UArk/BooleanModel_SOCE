from pysb.simulator import BngSimulator
import numpy as np
import matplotlib.pyplot as plt
import logging


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
            initials = [out_sp[-1] for out_sp in outputs[-1].species]
            # handle edge case where n_runs == 1 (note that n_runs is set to 1 below, so can't use 'if n_runs == 1')
            if len(np.array(initials).shape) == 1:
                initials = [outputs[-1].species[-1]]  # overwrite initials from above
            n_runs = 1  # number of sims based on number of initials after first step
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


class ObsToPlot(object):
    def __init__(self, name, label, color=None, marker=None):
        if isinstance(name, str):
            name = [name]
        self.name = name
        self.label = label
        self.color = color
        self.marker = marker


def plot_results(tspans, outputs, observables, multi_plots=False, save_plots=True, show_plots=False,
                 xlim=None, ylim=None, xlabel='iteration'):

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
            # make sure the observable exists. If not, skip it and send a warning
            if obs.name[0] not in out_observables.dtype.names:
                logging.warning("Observable '%s' not found. Skipping!" % obs.name[0])
                labels = labels[:-1]  # remove the last label since this observable doesn't exist
                break
            if len(np.array(out_observables[obs.name[0]]).shape) > 1:  # typical case -- n_runs > 1
                y = np.mean(out_observables[obs.name[0]], axis=0)
                if len(obs.name) > 1:
                    for i in range(1, len(obs.name)):
                        y += np.mean(out_observables[obs.name[i]], axis=0)
                    y /= len(obs.name)
            else:  # handle edge case where n_runs = 1
                y = np.array(out_observables[obs.name[0]])
                if len(obs.name) > 1:
                    for i in range(1, len(obs.name)):
                        y += np.array(out_observables[obs.name[i]])
                    y /= len(obs.name)
            line, = plt.plot(tspan, y, lw=3, color=obs.color)
            if step == 0:
                lines.append(line)
            # plot vertical lines for each step
            bottom, top = plt.ylim()
            plt.plot([tspan[-1], tspan[-1]], [bottom, top], 'k--')
    plt.xlabel(xlabel, fontsize=16)
    plt.ylabel('value', fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.legend(loc='upper left', handles=lines, labels=labels, bbox_to_anchor=(1, 1), fontsize=15)
    if xlim is not None:
        plt.xlim(xlim)
    if ylim is not None:
        plt.ylim(ylim)
    plt.tight_layout()
    if save_plots is not False:
        fname = "FIG_pysb_Boolean.pdf" if save_plots is True else save_plots
        plt.savefig(fname, format='pdf')
    if show_plots:
        plt.show()


if __name__ == '__main__':
    import os
    from pysb.importers.boolean import model_from_boolean

    mode = 'GSP'  # 'GSP', 'GA', 'ROA'
    n_runs = 100

    step_labels = [
        "equilibration",
        "remove external calcium and add BRAF inhibitor",
        "add pump inhibitor",
        "add external calcium"
    ]
    delta_ts = [50, 300, 10, 100]  # [50, 1000]  # [50, 50, 10, 100]
    conditions = [
        None,
        [("Ca_ext", False), ("BRAF_inhib", True)],
        [("Ca_pump_ER_inhib", True)],
        [("Ca_ext", True)]
    ]

    sim_steps = [SimStep(label, delta_t, condition) for label, delta_t, condition
                 in zip(step_labels, delta_ts, conditions)]

    for version in ['1_Pur']:

        print('mapk_soce_v%s.txt' % version)

        model = model_from_boolean(os.path.join('VERSIONS', 'mapk_soce_v%s.txt' % version), mode=mode)

        print('model created')

        tspans, outputs = sim_protocol(model, sim_steps, n_runs=n_runs, t_start=-sim_steps[0].delta_t, mode=mode,
                                       verbose=True)

        print('simulation complete')

        obs_names = [
            ['Ca_cyt_1_True_obs', 'Ca_cyt_2_True_obs', 'Ca_cyt_3_True_obs'],
            'BRAF_True_obs',
            'Ca_channel_True_obs',
            'Ca_ER_True_obs',
            'Ca_ext_True_obs',
            'Ca_pump_ER_True_obs',
            'ERK_True_obs',
            'Gene_Expr_True_obs',
            'MEK_True_obs'
        ]
        obs_labels = ['Ca_cyt', 'BRAF', 'Ca_channel', 'Ca_ER', 'Ca_ext', 'Ca_pump_ER', 'ERK', 'Gene_Expr', 'MEK']
        obs_colors = ['blue', 'green', 'black', 'red', 'purple', 'brown', 'yellow', 'orange', 'cyan']

        observables = [ObsToPlot(name, label, color) for name, label, color in zip(obs_names, obs_labels, obs_colors)]

        plot_results(tspans, outputs, observables, multi_plots=False,
                     save_plots='FIG_pysb_Boolean_v%s_8dayBRAFi.pdf' % version, show_plots=True,
                     xlim=(-5, 415), ylim=(-0.05, 1.05))
