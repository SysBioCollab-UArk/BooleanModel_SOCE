import boolean2
from boolean2 import util
import numpy as np
import pylab

n_delays = 50

initial_conditions = '''
BRAFi = False
BRAF = True
MEK = True
ERK = True
Gene_exp = False
Ca_ext = True
Ca_channel = True
Ca_cyt_1 = False
pumpi = False
Ca_pump_ER = True
Ca_ER = False
Ca_cyt_2 = False
Ca_cyt_3 = False
'''
# Add delay nodes
for i in range(n_delays):
    initial_conditions += 'Delay%d = False\n' % i

rules = '''
1: BRAF* = not BRAFi
1: MEK* = BRAF or (Ca_cyt_1 and Ca_cyt_2 and Ca_cyt_3)
1: ERK* = MEK
1: Ca_channel* = ERK or Gene_exp
1: Ca_cyt_1* = Ca_cyt_1 or (Ca_ext and Ca_channel) or Ca_ER
1: Ca_pump_ER* = not pumpi
1: Ca_ER* = Ca_cyt_1 and Ca_pump_ER
1: Ca_cyt_2* = Ca_cyt_1 and not Ca_pump_ER
1: Ca_cyt_3* = Ca_cyt_2 and Ca_ext and Ca_channel
'''
# Add delay rules
for i in range(n_delays):
    if i == 0:
        rules += '1: Delay0* = not ERK\n'
    else:
        rules += '1: Delay%d* = Delay%d\n' % (i, i-1)
    if i == n_delays-1:
        rules += '1: Gene_exp* = Delay%d\n' % (n_delays-1)

# Initial version of the model (for Step 1 below)
initial_model = initial_conditions + rules

# For plotting
species_to_plot = ['BRAF', 'Ca_channel', 'Ca_ER', 'Ca_ext', 'Ca_pump_ER', 'ERK', 'Gene_exp', 'MEK']
markers = ['o', '^', '*', 'd', 'H', 'v', '<', '>']
colors = ['green', 'black', 'red', 'purple', 'brown', 'yellow', 'orange', 'cyan']

# For storing trajectories
coll = util.Collector()

n_runs = 500  # number of Boolean runs

# Boolean update steps for each stage
equil_steps = 40  # 10
brafi_steps = 50  # 50
pumpi_steps = 10  # 10
ca_ext_steps = 100  # 10
pumpi_off_steps = 100  # 50

for i in range(n_runs):
    print(i)

    # Step 1, equilibration
    model = boolean2.Model(initial_model, mode='async')
    model.initialize()
    model.iterate(steps=equil_steps)
    coll.collect(states=model.states, nodes=model.nodes)

    # Step 2, add BRAF inhibitor
    initial_conditions = ""
    for species in model.states[-1].keys():
        if species == "BRAFi" or model.states[-1][species] is True:
            initial_conditions += '%s=True\n' % species
        else:
            initial_conditions += '%s=False\n' % species
    model = boolean2.Model(initial_conditions + rules, mode='async')
    model.initialize()
    model.iterate(steps=brafi_steps)
    coll_tmp = util.Collector()
    coll_tmp.collect(states=model.states, nodes=model.nodes)
    for species in coll.store.keys():
        coll.store[species][i] += coll_tmp.store[species][0][1:]

    # Step 3, remove external calcium AND add pump inhibitor
    initial_conditions = ""
    for species in model.states[-1].keys():
        if (species != "Ca_ext" and model.states[-1][species] is True) or species == "pumpi":
            initial_conditions += '%s=True\n' % species
        else:
            initial_conditions += '%s=False\n' % species
    model = boolean2.Model(initial_conditions + rules, mode='async')
    model.initialize()
    model.iterate(steps=pumpi_steps)
    coll_tmp = util.Collector()
    coll_tmp.collect(states=model.states, nodes=model.nodes)
    for species in coll.store.keys():
        coll.store[species][i] += coll_tmp.store[species][0][1:]

    # Step 4, add external calcium
    initial_conditions = ""
    for species in model.states[-1].keys():
        if species == "Ca_ext" or model.states[-1][species] is True:
            initial_conditions += '%s=True\n' % species
        else:
            initial_conditions += '%s=False\n' % species
    model = boolean2.Model(initial_conditions + rules, mode='async')
    model.initialize()
    model.iterate(steps=ca_ext_steps)
    coll_tmp = util.Collector()
    coll_tmp.collect(states=model.states, nodes=model.nodes)
    for species in coll.store.keys():
        coll.store[species][i] += coll_tmp.store[species][0][1:]

    # Step 5, remove pump inhibitor
    '''
    initial_conditions = ""
    for species in model.states[-1].keys():
        if species != "pumpi" and model.states[-1][species] is True:
            initial_conditions += '%s=True\n' % species
        else:
            initial_conditions += '%s=False\n' % species
    model = boolean2.Model(initial_conditions + rules, mode='async')
    model.initialize()
    model.iterate(steps=pumpi_off_steps)
    coll_tmp = util.Collector()
    coll_tmp.collect(states=model.states, nodes=model.nodes)
    for species in coll.store.keys():
        coll.store[species][i] += coll_tmp.store[species][0][1:]
    '''

# Get average node values
avgs = coll.get_averages(normalize=True)

# Calculate total cytosolic calcium levels
Ca_cyt_1 = np.array(avgs['Ca_cyt_1'])
Ca_cyt_2 = np.array(avgs['Ca_cyt_2'])
Ca_cyt_3 = np.array(avgs['Ca_cyt_3'])
Ca_cyt_avg = (Ca_cyt_1+Ca_cyt_2+Ca_cyt_3) / 3.0

# Plots
pylab.figure(figsize=(12.8, 4.8))
start = equil_steps
pylab.plot(Ca_cyt_avg[start:], 'sb-', label='Ca_cyt')
for sp, marker, color in zip(species_to_plot, markers, colors):
    pylab.plot(avgs[sp][start:], marker=marker, color=color, label=sp)
# pylab.legend(loc=0)
pylab.legend(loc="upper left", fontsize=15, ncol=1, bbox_to_anchor=(0.6, 0.9))  # (0.55, 0.8))
pylab.xlim(xmin=0)
pylab.ylim((-0.1, 1.1))
pylab.xlabel('iteration', fontsize=16)
pylab.ylabel('probability ON', fontsize=16)
pylab.xticks(fontsize=16)
pylab.yticks(fontsize=16)
pylab.tight_layout()

pylab.show()
