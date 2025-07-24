import pandas as pd
from plot_tune import *


# This is dummy EPIX code, replace with real file
# When ready real file should be a drop-in replacement
from EPIX_sample import get_EPICS_Tune, get_EPICS_Beta

def getValues() -> pd.DataFrame:

    plot_folder = 'Tune_Plots'
    cpymad_logfile = 'cpymad_logfile.txt'
    sequence_name = 'synchrotron'

    madx = cpymad_start(cpymad_logfile)

    # lattice_folder = '..\\Lattice_Files\\04_New_Harmonics\\'
    lattice_folder = '../Lattice_Files/04_New_Harmonics/'

    madx.call(file=lattice_folder+'ISIS.injected_beam')
    madx.call(file=lattice_folder+'ISIS.strength')
    madx.call(file=lattice_folder+'2023.strength')
    madx.call(file=lattice_folder+'ISIS.elements')
    madx.call(file=lattice_folder+'ISIS.sequence')

    cpymad_check_and_use_sequence(madx, cpymad_logfile, sequence_name)

    # Gets values for time
    time_periods = [-.4, -.2, 0, .5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8.1, 9, 9.8, 10.1]
    #time_periods = [0, 1]
    time_array = np.array(time_periods)

    # Gets Q values from EPIX
    qx_array = []
    qy_array = []

    for time in time_periods:
        EPIX_df = get_EPICS_Tune(time)
        qx_array.append(EPIX_df.loc[0, "Q_request"])
        qy_array.append(EPIX_df.loc[1, "Q_request"])

        
    # Check that arrays have same length
    assert(len(qx_array) == len(qy_array))
    assert(len(time_array) == len(qy_array))
    assert(len(time_array) == len(qx_array))

    # Calculating the currents
    tq_currents_df = tune_di_df(Qh=qx_array, Qv=qy_array, baseQh=4.331, baseQv=3.731, time_array=time_array, z=np.array([-4.73e-3, -5.99E-03, 4.45E-03, 2.40E-03]))
    #print(tq_currents_df)
    print('Expected values found')

    # Finally! Let's get the real tunes
    output_df = tune_calculation_iterator(madx, tq_currents_df, cpymad_logfile)
    #print(output_df)
    print('Actual values found')

    output = {
        "x":[],
        "y":[],
        "type":[],
        "time":[]
    }

    for i in range(0, len(output_df)):
        output['x'].append(output_df["Qh"][i])
        output['y'].append(output_df["Qv"][i])
        output['type'].append('set')
        output['time'].append(output_df["time"][i])

    for i in range(0, len(output_df)):
        output['x'].append(output_df["Machine Qh"][i])
        output['y'].append(output_df["Machine Qv"][i])
        output['type'].append('actual')
        output['time'].append(output_df["time"][i])

    df = pd.DataFrame(output)

    return df