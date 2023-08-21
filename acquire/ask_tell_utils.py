from dragonfly.opt.gp_bandit import CPGPBandit
import pandas as pd
import matplotlib.pyplot as plt
import circle_fit as cf
import numpy as np
from db_utils import add_electrolyte

def descriptive_ask(opt: CPGPBandit, all_considered_components):
    '''asks optimizer for a point, and returns dictionary with component formulas and amounts to test.

    Keyword arguments:
    all_considered_components -- same as all_chemicals from generate_config
    '''
    response = opt.ask()
    return {all_considered_components[i][0]:response[i] for i in range(len(all_considered_components))}
def csv_tell(ordered_component_dict:dict, opt: CPGPBandit, csv_filepath, x_cutoff, cell_constant = 10):
    '''inputs new electrolyte measurement result both into optimizer, and database.

    Keyword arguments:
    ordered_component_dict -- dictionary of formulas and amounts for electrolyte components; must be ordered in the correct way to tell cpgpbandit
    opt -- optimizer to tell result to
    csv_filepath -- the filepath to the eis nyquist csv data file for this run
    x_cutoff -- float, the x value where points to the right are not considered (used to remove warburg element)
    cell_constant -- float, cell constant for probe to derive conductivity from resistance.
    '''
    resistance, xc, yc = fit_arc_to_csv_data(csv_filepath,'/home/ubuntu/eis_db/acquire/data/plt.png',x_cutoff)
    conductivity = cell_constant/resistance
    opt.tell([(ordered_component_dict.values(), conductivity)])
    add_electrolyte(ordered_component_dict,conductivity,-1,-1)

def fit_arc_to_csv_data(csv_path, output_filepath, x_cutoff: float):
    '''fits arc to csv data and returns diameter, and x,y coordinates for center.

    Keyword arguments:
    csv_path -- filepath to csv nyquist data
    output_filepath -- if graph code is uncommented, output filepath for graph image.
    x_cutoff -- float, the x value where points on the right are not considered (used to remove warburg element)
    '''
    df = pd.read_csv(csv_path)
    x_data = df['ZRe/Ohm'].values
    y_data = df['negZIm/Ohm'].values
    
    #filter data based on x_cutoff
    if x_cutoff is not None:
        mask = x_data <= x_cutoff
        x_data = x_data[mask]
        y_data = y_data[mask]

    data_combined = np.vstack((x_data, y_data)).T

    xc, yc, radius, _ = cf.least_squares_circle(data_combined)
    '''
    plt.scatter(x_data, y_data, label='Data', color='blue')

    theta = np.linspace(0, 2*np.pi, 100)
    x_fit = xc + radius * np.cos(theta)
    y_fit = yc + radius * np.sin(theta)
    plt.plot(x_fit, y_fit, label='Fitted Arc', color='red')
    
    plt.xlabel('ZRe/Ohm')
    plt.ylabel('negZIm/Ohm')
    plt.legend()
    plt.title('Arc Fitted on Noisy Data')
    plt.grid(True)
    plt.gca().set_aspect('equal', adjustable='box')
    plt.savefig(output_filepath)
    plt.show()
    '''

    return 2*radius, xc, yc

filepath = '/home/ubuntu/eis_db/acquire/data/eis_circle_fit_test_data.csv'
img_filepath='/home/ubuntu/eis_db/acquire/data/plt.png'
fit_arc_to_csv_data(filepath,img_filepath,65.85)