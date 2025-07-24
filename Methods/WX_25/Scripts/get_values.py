from plot_tune import *
import pandas as pd

def getValues():

    plot_folder = 'Tune_Plots'
    cpymad_logfile = 'cpymad_logfile.txt'
    sequence_name = 'synchrotron'

    madx = cpymad_start(cpymad_logfile)

    lattice_folder = 'ISIS_Synchrotron_Model\\Lattice_Files\\04_New_Harmonics\\'

    madx.call(file=lattice_folder+'ISIS.injected_beam')
    madx.call(file=lattice_folder+'ISIS.strength')
    madx.call(file=lattice_folder+'2023.strength')
    madx.call(file=lattice_folder+'ISIS.elements')
    madx.call(file=lattice_folder+'ISIS.sequence')

    cpymad_check_and_use_sequence(madx, cpymad_logfile, sequence_name)

    # dummy values for set tunes (get values from EPICS later)
    qx_array = [4.315, 4.270, 4.270, 4.250, 4.235, 4.205, 4.170, 4.190, 4.18, 4.18, 4.18, 4.17, 4.165, 4.165, 4.165, 4.18, 4.18, 4.175]
    qy_array = [3.82, 3.82, 3.81, 3.805, 3.800, 3.825, 3.680, 3.680, 3.69, 3.7, 3.7, 3.695, 3.695, 3.695, 3.692, 3.69, 3.685, 3.665]

    # dummy values for time array
    time_array = np.array([-0.1, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0, 10.0])


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


    ## Plotly
    
    import seaborn as sns
    import plotly.express as px

    # plotting set & actual tunes
    fig = px.scatter(df,
        x="x", 
        y="y", 
        color="time", 
        symbol="type",
    )
    fig.update_layout(
        title="Set & Actual Tunes",
        xaxis_title='Qh',
        yaxis_title='Qv',
        legend=dict(x=0, y=1, traceorder='normal', orientation='h')
    )
    fig.show()

    return df

print(getValues())
# output_df = tune_calculation_iterator(madx, tq_currents_df, cpymad_logfile)


