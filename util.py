from pysb_Boolean import SimStep
import re
from io import StringIO

base_mapk_soce_model = '''
BRAF = True
MEK = True
ERK = True
Ca_ext = True
Ca_channel = True
Ca_pump_ER = True
Ca_ER = True
Ca_cyt_1 = True
Ca_cyt_2 = False
Ca_cyt_3 = False

1: MEK* = BRAF
1: ERK* = MEK
1: Ca_channel* = ERK
1: Ca_ER* = (Ca_cyt_1 or Ca_cyt_2 or Ca_cyt_3) and Ca_pump_ER
1: Ca_cyt_1* = (Ca_cyt_1 or (Ca_ext and Ca_channel) or Ca_ER) and not Ca_cyt_2 and not Ca_cyt_3 or (Ca_cyt_2 and not Ca_cyt_3 or Ca_cyt_3 and not Ca_cyt_2) and not Ca_pump_ER or (Ca_cyt_2 and Ca_cyt_3 and Ca_ext and Ca_channel)
1: Ca_cyt_2* = (Ca_cyt_2 or (Ca_ext and Ca_channel) or Ca_ER) and not Ca_cyt_1 and not Ca_cyt_3 or (Ca_cyt_1 and not Ca_cyt_3 or Ca_cyt_3 and not Ca_cyt_1) and not Ca_pump_ER or (Ca_cyt_1 and Ca_cyt_3 and Ca_ext and Ca_channel)
1: Ca_cyt_3* = (Ca_cyt_3 or (Ca_ext and Ca_channel) or Ca_ER) and not Ca_cyt_1 and not Ca_cyt_2 or (Ca_cyt_1 and not Ca_cyt_2 or Ca_cyt_2 and not Ca_cyt_1) and not Ca_pump_ER or (Ca_cyt_1 and Ca_cyt_2 and Ca_ext and Ca_channel)
'''


def create_mapk_soce_model(base_model, new_nodes=None, new_edges=None, inhibited_nodes=None):

    # extract initial conditions and rules from the base model
    initial_conditions = ''
    rules = ''
    f = StringIO(base_model)
    for line in f.readlines():
        if re.search(r'^\w+\s*=\s*\w+\n', line):
            initial_conditions += line
        elif re.search(r'^\d+:\s*\w+\*\s*=\s*.*\n', line):
            rules += line

    if new_nodes is None:
        new_nodes = []
    if new_edges is None:
        new_edges = []
    if inhibited_nodes is None:
        inhibited_nodes = []

    # modify initial conditions
    for node in new_nodes:
        initial_conditions += '%s = %s\n' % (node['name'], node['initial'])

    n = 0
    for edge in new_edges:
        if edge['n_delays'] > 0:
            for j in range(edge['n_delays']):
                initial_conditions += 'Delay%d_%d = False\n' % (n, j)
            n += 1

    for node in inhibited_nodes:
        initial_conditions += '%s_inhib = False\n' % node

    # modify rules
    for node in new_nodes:
        if node['rule'] is not None:
            rules += '1: %s* = %s\n' % (node['name'], node['rule'])

    n = 0
    for edge in new_edges:
        final_condition = edge['condition']  # this lets us create the correct rule if there are zero delays
        # create all delay rules
        if edge['n_delays'] > 0:
            for j in range(edge['n_delays']):
                delay_prefix = 'Delay%d_' % n
                if j == 0:
                    rules += '1: %s0* = %s\n' % (delay_prefix, edge['condition'])
                else:
                    rules += '1: %s%d* = %s%d\n' % (delay_prefix, j, delay_prefix, j - 1)
                final_condition = '%s%d' % (delay_prefix, j)
            n += 1
        # create last rule in the chain
        target_pattern = r'(%s\*\s*=\s*)(.+)' % edge['target']
        m = re.search(target_pattern, rules)
        if m is not None:  # a rule for this node already exists
            rules = re.sub(target_pattern, r'%s or %s' % (m.group(0), final_condition), rules)
        else:
            rules += '1: %s* = %s\n' % (edge['target'], final_condition)

    for node in inhibited_nodes:
        pattern = r'(%s\*\s*=\s*)(.+)' % node
        m = re.search(pattern, rules)
        if m is None:
            rules += '1: %s* = not %s_inhib\n' % (node, node)
        else:
            rules = re.sub(pattern, r'%s(%s) and not %s_inhib' % (m.group(1), m.group(2), node), rules)

    return initial_conditions + rules


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
            [("Ca_pump_ER_inhib", True)],
            [("Ca_ext", True)]
        ]

    # 30 min BRAFi
    elif which == "30_min_BRAFi":

        step_labels = [
            "equilibration",
            "remove external calcium and add BRAF inhibitor",
            "add pump inhibitor",
            "add external calcium"
        ]
        delta_ts = [50, 50, 10, 100]
        conditions = [
            None,
            [("Ca_ext", False), ("BRAF_inhib", True)],
            [("Ca_pump_ER_inhib", True)],
            [("Ca_ext", True)]
        ]

    # 8 days BRAFi
    elif which == "8_days_BRAFi":

        step_labels = [
            "equilibration",
            "remove external calcium and add BRAF inhibitor",
            "add pump inhibitor",
            "add external calcium"
        ]
        delta_ts = [50, 300, 10, 100]
        conditions = [
            None,
            [("Ca_ext", False), ("BRAF_inhib", True)],
            [("Ca_pump_ER_inhib", True)],
            [("Ca_ext", True)]
        ]

    # 8 days BRAFi + 15 min MEKi
    elif which == "8_days_BRAFi_plus_MEKi":

        step_labels = [
            "equilibration",
            "remove external calcium and add BRAF inhibitor",
            "add MEK inhibitor",
            "add pump inhibitor",
            "add external calcium"
        ]
        delta_ts = [50, 300, 10, 10, 100]
        conditions = [
            None,
            [("Ca_ext", False), ("BRAF_inhib", True)],
            [("MEK_inhib", True)],
            [("Ca_pump_ER_inhib", True)],
            [("Ca_ext", True)]
        ]

    return [SimStep(label, dt, cond) for label, dt, cond in zip(step_labels, delta_ts, conditions)]


if __name__ == '__main__':
    import os
    model_versions = list()

    # version 1
    model_versions.append(
        create_mapk_soce_model(
            base_mapk_soce_model,
            new_nodes=[{'name': 'Gene_Expr', 'initial': False, 'rule': None}],
            new_edges=[{'target': 'Ca_channel', 'condition': 'Gene_Expr', 'n_delays': 0},
                       {'target': 'MEK', 'condition': '(Ca_cyt_1 and Ca_cyt_2 and Ca_cyt_3)', 'n_delays': 0},
                       {'target': 'Gene_Expr', 'condition': 'not ERK', 'n_delays': 100}],
            inhibited_nodes=['BRAF', 'MEK', 'Ca_pump_ER']
        )
    )

    # version 2
    model_versions.append(
        create_mapk_soce_model(
            base_mapk_soce_model,
            new_edges=[{'target': 'MEK',
                        'condition': '(Ca_cyt_1 and Ca_cyt_2) or (Ca_cyt_1 and Ca_cyt_3) or (Ca_cyt_2 and Ca_cyt_3)',
                        'n_delays': 100}],
            inhibited_nodes=['BRAF', 'MEK', 'Ca_pump_ER']
        )
    )

    # version 3
    model_versions.append(
        create_mapk_soce_model(
            base_mapk_soce_model,
            new_nodes=[{'name': 'Gene_Expr', 'initial': False, 'rule': None}],
            new_edges=[{'target': 'Ca_channel', 'condition': 'Gene_Expr', 'n_delays': 0},
                       {'target': 'MEK', 'condition': '(Ca_cyt_1 and Ca_cyt_2 and Ca_cyt_3)', 'n_delays': 0},
                       {'target': 'Gene_Expr', 'condition': 'not ERK', 'n_delays': 100},
                       {'target': 'Ca_pump_ER', 'condition': 'MEK', 'n_delays': 100}],
            inhibited_nodes=['BRAF', 'MEK', 'Ca_pump_ER']
        )
    )

    # version 4
    model_versions.append(
        create_mapk_soce_model(
            base_mapk_soce_model,
            new_edges=[{'target': 'MEK',
                        'condition': '(Ca_cyt_1 and Ca_cyt_2) or (Ca_cyt_1 and Ca_cyt_3) or (Ca_cyt_2 and Ca_cyt_3)',
                        'n_delays': 100},
                       {'target': 'Ca_pump_ER', 'condition': 'MEK', 'n_delays': 100}],
            inhibited_nodes=['BRAF', 'MEK', 'Ca_pump_ER']
        )
    )

    # version 5
    model_versions.append(
        create_mapk_soce_model(
            base_mapk_soce_model,
            new_nodes=[{'name': 'Gene_Expr', 'initial': False, 'rule': None}],
            new_edges=[{'target': 'Ca_channel', 'condition': 'Gene_Expr', 'n_delays': 0},
                       {'target': 'MEK',
                        'condition': '(Ca_cyt_1 and Ca_cyt_2) or (Ca_cyt_1 and Ca_cyt_3) or (Ca_cyt_2 and Ca_cyt_3)',
                        'n_delays': 0},
                       {'target': 'Gene_Expr', 'condition': 'not ERK', 'n_delays': 100},
                       {'target': 'Ca_pump_ER', 'condition': 'MEK', 'n_delays': 100}],
            inhibited_nodes=['BRAF', 'MEK', 'Ca_pump_ER']
        )
    )

    # version 6
    model_versions.append(
        create_mapk_soce_model(
            base_mapk_soce_model,
            new_nodes=[{'name': 'Gene_Expr', 'initial': False, 'rule': None}],
            new_edges=[{'target': 'Ca_channel', 'condition': 'Gene_Expr', 'n_delays': 0},
                       {'target': 'MEK',
                        'condition': '(Ca_cyt_1 and Ca_cyt_2) or (Ca_cyt_1 and Ca_cyt_3) or (Ca_cyt_2 and Ca_cyt_3)',
                        'n_delays': 100},
                       {'target': 'Gene_Expr', 'condition': 'not ERK', 'n_delays': 100},
                       {'target': 'Ca_pump_ER', 'condition': 'MEK', 'n_delays': 100}],
            inhibited_nodes=['BRAF', 'MEK', 'Ca_pump_ER']
        )
    )

    for i, version in enumerate(model_versions):
        with open(os.path.join('VERSIONS', 'mapk_soce_v%d.txt' % (i + 1)), "w") as outfile:
            outfile.write(version)
