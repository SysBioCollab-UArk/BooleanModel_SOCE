# importing libraries
import pylab  # to plot
import os  # to direct path
import boolean2  # for simulations
from boolean2 import Model
from boolean2 import util  # to collect data from async boolean simulations

# Defining the path to the folder to save plots and data

main_folder = "C:\Users\harsi\Desktop\BRAFi_SOCE figures"
path_to_save = main_folder + "/figure_9"

# if the path does not exist this code will generate a folder
isExist = os.path.exists(path_to_save)
if not isExist:
    os.mkdir(path_to_save)

# text for initial condition of each node and the interaction between
# node in the model.

text = """
BRAFi = True
BRAF = False
MEK = False
ERK = False
Ca_channel = True
Gen_exp = False
ca_cyt_normal= True
Ca_ER = True
Ca_ext = True
Ca_er_pump = False
Ca_er_pumpi = True
Ca_cyt_high = False

1: BRAF* = not BRAFi
1: MEK* = BRAF
1: ERK* = MEK or Ca_cyt_high
1: Ca_channel* = ERK or Gen_exp
1: Gen_exp* = not ERK
1: ca_cyt_normal = Ca_ext and Ca_channel or Ca_ER
1: Ca_cyt_high* = ca_cyt_normal and not Ca_ER
1: Ca_ER* = ca_cyt_normal and Ca_er_pump
1: Ca_er_pump* = not Ca_er_pumpi

"""

# if mode is sync uncomment this part

# model = Model(text=text, mode='sync')
# model.initialize()
# model.iterate(steps=10)
#
# # for node in model.data:
# #     print(node, model.data[node])
#
# p1 = pylab.plot(model.data["Ca_cyt"], 'sb-', label='Ca_cyt')
# p2 = pylab.plot(model.data["ERK"], 'vr-', label='ERK')
# p3 = pylab.plot(model.data["BRAFi"], 'og-', label='BRAFi')
# pylab.legend()
# pylab.ylim((-0.1, 1.1))
# pylab.show()

# if mode is async uncomment this part

# boolean simulation
coll = util.Collector()
for i in range(1000):
    model = boolean2.Model(text, mode='async')
    model.initialize()
    model.iterate(steps=10)

    # takes all nodes
    nodes = model.nodes
    coll.collect(states=model.states, nodes=nodes)

# Taking average of the states for 1000 runs

avgs = coll.get_averages(normalize=True)

# plotting the data

p1 = pylab.plot(avgs['Ca_cyt'], 'sb-', label='Ca_cyt')  # type: object
p2 = pylab.plot(avgs['BRAF'], 'vr-', label='BRAF')
p3 = pylab.plot(avgs['ERK'], 'og-', label='ERK')
p4 = pylab.plot(avgs['Ca_cyt_high'], 'sy-', label='Ca_cyt_high')
pylab.legend()
pylab.ylim((-0.1, 1.1))
pylab.savefig(path_to_save + '/long-term_BRAF_treatment_ER_pump_off.png')
pylab.show()

# write a function to save the figures
# text_file = open(path_to_save+'/Param.txt', 'w')
# n = text_file.write(text)
# text_file.close()
