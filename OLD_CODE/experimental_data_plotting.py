import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal

path_to_folder = "C:/Users/harsi/Downloads/08_2021_data/08_2021_data"


'''
This function reads the excel data files sheet Dara (ROIs)
'''
def reading_data_files(path_to_file):
    data = pd.read_excel(path_to_file, sheet_name = 'Data (ROIs)' )

    return data

'''
Mean and standard deviation for each time point for all the cells.
'''
def estimate_mean_sdev(data):
    stats_dat = pd.DataFrame()

    stats_dat['mean'] = data.iloc[:, 5:].mean(axis=1)
    stats_dat['std_dev'] = data.iloc[:, 5:].std(axis=1)

    data['Time [s]'] = data['Time [s]']/60
    stats_dat = stats_dat.set_index(data['Time [s]'])

    return stats_dat

'''
This function plots the mean and standard deviation for each experiment.
It also estimates the local maxima's stored in pos variable.
These local maxima's are marked in red on plots
and are used for bar plots in next function.
'''
def line_plot(stats, path_to_folder, exp):
    x = np.array(stats.index)
    y = np.array(stats['mean'])
    ci = np.array(stats['std_dev'])

    pos = signal.argrelextrema(y, np.greater, order=25) # estimating the local maxima
    dat = stats.iloc[pos]
    dat_bar = pd.DataFrame()

    if np.size(pos)>3 :
        dat_bar = dat.iloc[1:-1]

    elif np.size(pos)>2 :
        dat_bar = dat.iloc[1:]

    print(dat_bar)

    fig, ax = plt.subplots()
    ax.plot(x, y)
    ax.fill_between(x, (y - ci), (y + ci), color='b', alpha=.1)
    ax.scatter(x[pos], dat['mean'], color='r')
    ax.set_xlim(0,27.5)
    ax.set_ylim(0.5,1.2)

    plt.savefig(path_to_folder + '/' + exp + '.png')
    plt.close()
    return dat_bar

def plotting_bar_plot(data_mean, data_sdev,path_to_folder):

    fig = plt.figure()
    ax = plt.subplots(111)
    ax = data_mean.plot(kind = 'bar', yerr = data_sdev)
    plt.legend(loc=2, prop={'size': 6})

    # plt.show()

    plt.savefig(path_to_folder + '/bar_plot.png', dpi = 400)



categories = ['naive', 'plx24hr_constant', 'plx24hours_holiday', 'plx3day_holiday' ,
              'plx3day_constant', 'plx8day_constant', 'plx8day_holiday']

data_bar_plot_mean = pd.DataFrame()
data_bar_plot_sdev = pd.DataFrame()

for i in categories:
    exp = 'A375_' + i
    filename = exp + '.xlsx'
    path_to_file = path_to_folder + '/' + filename
    data = reading_data_files(path_to_file)
    stats = estimate_mean_sdev(data)
    local_maxima = line_plot(stats, path_to_folder, exp)
    local_maxima = local_maxima.reset_index()
    data_bar_plot_mean[exp] = local_maxima['mean']
    data_bar_plot_sdev[exp] = local_maxima['std_dev']


plotting_bar_plot(data_bar_plot_mean,data_bar_plot_sdev,path_to_folder)