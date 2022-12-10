from pysb.importers.boolean import model_from_boolean
from pysb.simulator import BngSimulator
import numpy as np
import matplotlib.pyplot as plt

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


observables= [
    'BRAF_True_obs',
    'ERK_True_obs',
    'Ca_cyt_1_True_obs',
    'Ca_cyt_2_True_obs',
    'Ca_cyt_3_True_obs',
    'Ca_channel_True_obs',
    'Ca_ER_True_obs',
    'Ca_pump_ER_True_obs']


tspan=np.linspace(0,10,101)
sim=BngSimulator(model,tspan,verbose=False)



n_runs=10
traj=[None]*n_runs


for n in range(n_runs):
    print(n)
    #Step 1 equilibration
    output=sim.run()
    # for sp in model.species:
    #     print(sp)
    # quit()
    traj[n]= output.observables['BRAF_True_obs']

    #step2 Adding Brafi
    initials=output.species[-1]
    initials[0]=0
    initials[1]=1
    output2=sim.run(initials=initials)
    traj[n]=np.append(traj[n],output2.observables['BRAF_True_obs'])


    # Step 3, remove external calcium AND add pump inhibitor

    # Step 4, add external calcium

    # Step 5, remove pump inhibitor


avg_traj=np.mean(np.array(traj),axis=0)
#avg_traj=np.mean(np.array([output.observables[i]['Ca_channel_True_obs'] for i in range(n_runs)]),axis=0)
time=np.append(tspan,tspan+tspan[-1])
plt.plot(time,avg_traj,lw=2,label='BRAF_True_obs')




plt.legend(loc=0)
plt.show()
quit()


# BRAF,ERK,Ca_cyt,Ca_channel,Ca_ER,Ca_Punmp_ER
observables= [
    'BRAF_True_obs',
    'ERK_True_obs',
    'Ca_cyt_1_True_obs',
    'Ca_cyt_2_True_obs',
    'Ca_cyt_3_True_obs',
    'Ca_channel_True_obs',
    'Ca_ER_True_obs',
    'Ca_pump_ER_True_obs']

for obs_name in observables:
    plt.plot(tspan,output.observables[obs_name],lw=2,label=obs_name)

plt.legend(loc=0)
plt.show()
