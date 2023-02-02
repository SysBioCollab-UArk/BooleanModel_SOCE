from pysb.importers.boolean import model_from_boolean
from pysb.simulator import BngSimulator
import numpy as np
import matplotlib.pyplot as plt
from pysb.bng import *
from pysb.generator.bng import BngPrinter
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
n_runs=100

observables= [
    ('BRAF_True_obs', 'b'),
    ('ERK_True_obs','r'),
    ('Ca_cyt_1_True_obs','g'),
    ('Ca_cyt_2_True_obs','k'),
    ('Ca_cyt_3_True_obs','c'),
    ('Ca_channel_True_obs','m'),
    ('Ca_ER_True_obs','y'),
    ('Ca_pump_ER_True_obs','0.5'),
    ('MEK_True_obs','purple'),
    ('Gene_exp_True_obs','orange')
]

traj={}

# for obs in observables:
#     traj[obs]=[None]*n_runs


#tspan1 = np.linspace(0, 10, 101)
# sim=BngSimulator(model,verbose=True)
# sim.run(tspan=tspan1)
# quit()

nsims=10
with BngFileInterface(model) as bng:
    bng.action('generate_network')
    for n in range(nsims):
        print (n)
        # step1 = equilibration
        tspan1 = np.linspace(0, 10, 101)
        tspan_all = tspan1
        bng.action('simulate', method='ssa', t_start=tspan_all[0],t_end=tspan_all[-1] ,n_steps=len(tspan1)-1)

        #step2= adding BRAFi
        tspan2 = np.linspace(0, 10, 101)
        bng.set_concentration("BRAFi(state~False)",0)
        bng.set_concentration("BRAFi(state~True)",1)
        kwarg={'continue': 1}
        bng.action('simulate', method='ssa',t_start=tspan_all[-1],t_end=tspan_all[-1]+tspan2[-1] ,n_steps=len(tspan2)-1, **kwarg)
        tspan_all = np.append(tspan_all, tspan2[1:] + tspan1[-1])

        #step3= remove external calcium and adding pump inhibitor
        tspan3 = np.linspace(0, 10, 101)
        bng.set_concentration("Ca_ext(state~False)",1)
        bng.set_concentration("Ca_ext(state~True)",0)
        bng.set_concentration("pumpi(state~False)",0)
        bng.set_concentration("pumpi(state~True)",1)
        bng.action('simulate',method='ssa',t_start=tspan_all[-1],t_end=tspan_all[-1]+tspan3[-1],n_steps=len(tspan3)-1,**kwarg)
        tspan_all = np.append(tspan_all, tspan3[1:] + tspan1[-1] + tspan2[-1])

        #step4=add external calcium
        tspan4 = np.linspace(0, 10, 101)
        bng.set_concentration("Ca_ext(state~False)",0)
        bng.set_concentration("Ca_ext(state~True)",1)
        bng.action('simulate',method='ssa',t_start=tspan_all[-1],t_end=tspan_all[-1]+tspan4[-1],n_steps=len(tspan4)-1,**kwarg)
        tspan_all = np.append(tspan_all, tspan4[1:] + tspan1[-1] + tspan2[-1] + tspan3[-1])

        #step5=
        # remove pump inhibitor
        tspan5 = np.linspace(0, 10, 101)
        bng.set_concentration("pumpi(state~False)",1)
        bng.set_concentration("pumpi(state~True)",0)
        bng.action('simulate',method='ssa',t_start=tspan_all[-1],t_end=tspan_all[-1]+tspan5[-1],n_steps=len(tspan5)-1,**kwarg)
        tspan_all = np.append(tspan_all, tspan5[1:] + tspan1[-1] + tspan2[-1] + tspan3[-1] + tspan4[-1])

        #reset concentration
        bng.action('resetConcentrations')

    #run the simulation and plot the results
    bng.execute()
    result=bng.read_simulation_results()
    print (result.shape)
    for obs in observables:
        plt.figure(obs[0])
        plt.plot(tspan_all, result[obs[0]],lw=2,color=obs[1],label=obs[0])


plt.show()
quit()

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
