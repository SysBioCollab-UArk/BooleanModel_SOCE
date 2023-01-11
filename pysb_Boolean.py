from pysb.importers.boolean import model_from_boolean
from pysb.simulator import BngSimulator
import numpy as np
import matplotlib.pyplot as plt
from pysb.bng import BngFileInterface

model = model_from_boolean('mapk_soce.txt', mode='GSP')

# print(model)
# print (model.monomers)
# print (model.rules)
# print (model.parameters)
# print (len(model.parameters_rules()))
# print(model.observables)

# x=[[1,2,3],
#    [4,5,6],
#    [7,8,9]]
#
# print(np.mean(x,axis=0))
# quit()
n_runs=10

observables= [
    'BRAF_True_obs',
    'ERK_True_obs',
    'Ca_cyt_1_True_obs',
    'Ca_cyt_2_True_obs',
    'Ca_cyt_3_True_obs',
    'Ca_channel_True_obs',
    'Ca_ER_True_obs',
    'Ca_pump_ER_True_obs',
    'MEK_True_obs',
    'Gene_exp_True_obs'
]

traj={}

for obs in observables:
    traj[obs]=[None]*n_runs

# sim = BngSimulator(model, verbose=True)

with BngFileInterface(model) as bng:
    bng.action('generate_network')

    for n in range(n_runs):
        print(n)
        # Step 1 equilibration
        # tspan1 = np.linspace(0, 10, 101)
        # output = sim.run(tspan=tspan1)
        bng.action('simulate', method='ssa', t_end=10, n_steps=100)
        bng.execute()
        output = bng.read_simulation_results()
        # TODO: still need to figure out how to use the BngFileInterface correctly
        # for sp in model.species:
        #     print(sp)
        # quit()
        for obs in observables:
            traj[obs][n]=output.observables[obs] # here we need the nth simulation result for each observable

        #step2 Adding Brafi
        tspan2=np.linspace(0, 10, 101)
        initials=output.species[-1]
        initials[0]=0
        initials[1]=1
        output2 = sim.run(tspan=tspan2,initials=initials)
        for obs in observables:
            traj[obs][n]=np.append(traj[obs][n],output2.observables[obs])


        # Step 3, remove external calcium AND add pump inhibitor


        # Step 4, add external calcium

        # Step 5, remove pump inhibitor

avg_traj={}
for obs in observables:
    avg_traj[obs]=np.mean(np.array(traj[obs]),axis=0)

    time=np.append(tspan1,tspan2+tspan1[-1])
    plt.plot(time,avg_traj[obs],lw=2,label=obs)
plt.legend(loc=0)
plt.show()
quit()

# BRAF,ERK,Ca_cyt,Ca_channel,Ca_ER,Ca_Punmp_ER
for obs_name in observables:
    plt.plot(tspan,output.observables[obs_name],lw=2,label=obs_name)

plt.legend(loc=0)
plt.show()
