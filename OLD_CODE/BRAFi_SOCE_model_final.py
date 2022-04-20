'''
This pipeline requires model_rules_0.txt file with
the rules and initial conditions to run.
'''

# importing all the essential libraries

import matplotlib.pylab as pylab
import os
import boolean2
import pandas as pd
from boolean2 import util

# Creating a path to save plots, results and the parameters

main_folder = os.path.join(os.getcwd(), "BRAFi_SOCE_figures")
if not os.path.exists(main_folder):
    os.mkdir(main_folder)

path_to_save = os.path.join(main_folder, "Network_3")
if not os.path.exists(path_to_save):
    os.mkdir(path_to_save)

# function to simulate model

def async_model_sim(model_rules, iterations):
    coll = util.Collector()
    for i in range(1000):
        model = boolean2.Model(model_rules, mode='async')
        model.initialize()
        model.iterate(steps=iterations)

        # takes all nodes
        nodes = model.nodes
        coll.collect(states=model.states, nodes=nodes)

    avgs = coll.get_averages(normalize=True)
    df = pd.DataFrame.from_dict(avgs, orient="index")
    # df.to_csv(path_to_save)

    return (avgs, df, nodes)


'''
This function updates the initial condition for the next set of simulations
with the end point of the states in the last simulation. A new parameter file is
created and saved in the folder path defined above.
'''

def updating_initial_condition(rules, nodes, df, path_to_save, count):
    with open("model_rules_" + str(count) + ".txt", 'w') as f:
        f.write(rules)

    with open("model_rules_" + str(count) + ".txt", 'r') as f:
        lines = f.readlines()

    for node, i in zip(nodes, range(0, 23)):
        x = df.loc[node, 10]
        # print(x)
        x = int(x)
        n_line = str(node) + ' = ' + str(bool(x)) + '\n'
        lines[i] = n_line
        # # print(lines[i])

    with open("model_rules_" + str(count) + ".txt", 'w') as f:
        f.writelines(lines)

'''
This function generates the plots for the simulations. 
'''
def plots_for_simulation(average_states, count, path_to_save):
    ca_levels = [sum(x) for x in zip(average_states['Ca_cyt_normal'], average_states['Ca_ER_dump'],
                                     average_states['Ca_SOCE'])]
    ca_levels = [x / 3 for x in ca_levels]

    p1 = pylab.plot(ca_levels, 'sb-', label='Cytoplasmic Ca')  # type: object
    p2 = pylab.plot(average_states['BRAF'], 'vr-', label='BRAF')
    p3 = pylab.plot(average_states['ERK'], 'og-', label='ERK')
    p4 = pylab.plot(average_states['MEK'], 'vc-', label='MEK')
    p5 = pylab.plot(average_states['UK_node'], 'sy-', label='UK_node')
    pylab.legend()
    pylab.ylim((-0.1, 1.1))
    pylab.savefig(os.path.join(path_to_save, 'short_fig_' + str(count) + '.png'))
    pylab.close()


'''
Following lines of code run the simulations for different experimental condition.
Final states for all the runs for different experimental condition are saved in combined_states dictionary.
boolean2.modify_states updates the nodes that stay off/on for that set of simulations.
Depending upon the treatment duration short/long the steps for BRAFi treatment are determined.
ERK_state is the mean of ERK state for the last 20 steps is determined for the next initial condition run. 
Only if it is 0, the UK_node is turned on otherwise it stays off for the rest of simulations. 
'''
count = range(5)
treatment = 'short'
combined_states = {}
for i in count:
    with open("model_rules_" + str(i) + ".txt", 'r') as f:
        model_rules = f.read()

    if (i < 1):

        on = []
        off = ['BRAFi', 'Ca_ER_pumpi ']
        new_txt = boolean2.modify_states(model_rules, turnon=on, turnoff=off)
        average_states, df, nodes = async_model_sim(new_txt, iterations=10)
        updating_initial_condition(model_rules, nodes, df, path_to_save, count=i+1)
        plots_for_simulation(average_states,i,path_to_save)
        combined_states = average_states



    elif(i==1):
        if treatment == 'short':
            steps = 10
        else:
            steps = 30
        on = ['BRAFi']
        off = []
        new_txt = boolean2.modify_states(model_rules, turnon=on, turnoff=off)
        average_states, df, nodes = async_model_sim(new_txt,iterations=steps)
        updating_initial_condition(model_rules, nodes, df, path_to_save, count=i+1)
        plots_for_simulation(average_states,i,path_to_save)
        for key, value in combined_states.iteritems():
            value.extend(average_states[key])

    elif(i == 2):
        ERK_state = combined_states['ERK']
        ERK_state_mean = (sum(ERK_state[-20:])/len(ERK_state[-20:]))
        if ERK_state_mean == 0:
            on = ['UK_node']
        off = ['Ca_ext']
        new_txt = boolean2.modify_states(model_rules, turnon=on, turnoff=off)
        average_states, df, nodes = async_model_sim(new_txt, iterations=10)
        updating_initial_condition(model_rules, nodes, df, path_to_save, count=i + 1)
        plots_for_simulation(average_states,i,path_to_save)
        for key, value in combined_states.iteritems():
            value.extend(average_states[key])

    elif(i == 3):
        on = ['Ca_ER_pumpi ']
        off = []
        new_txt = boolean2.modify_states(model_rules, turnon=on, turnoff=off)
        average_states, df, nodes = async_model_sim(new_txt,iterations=10)
        updating_initial_condition(model_rules, nodes, df, path_to_save, count=i+1)
        plots_for_simulation(average_states,i,path_to_save)
        for key, value in combined_states.iteritems():
            value.extend(average_states[key])

    else:
        on = ['Ca_ext']
        off = []
        new_txt = boolean2.modify_states(model_rules, turnon=on, turnoff=off)
        average_states, df, nodes = async_model_sim(new_txt, iterations=30)
        updating_initial_condition(model_rules, nodes, df, path_to_save, count=i + 1)
        # plots_for_simulation(average_states,i,path_to_save)
        for key, value in combined_states.iteritems():
            value.extend(average_states[key])

# print(combined_states)

plots_for_simulation(combined_states, "combined", path_to_save)
# plots_for_simulation_delay(combined_states, "combined_2", path_to_save)
