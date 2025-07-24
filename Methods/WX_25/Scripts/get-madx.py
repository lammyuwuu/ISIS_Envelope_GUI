from plot_tune import *
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

max_E = 800 # 800 MeV
cycle_time = 0.0 # 1.5 milliseconds into the 10 ms acceleration cycle of the synchrotron
cpymad_set_isis_cycle_time(madx, max_E, cycle_time)

twiss_0 = cpymad_madx_twiss(madx, cpymad_logfile, sequence_name)

print(twiss_0)







# output_df = tune_calculation_iterator(madx, tq_currents_df, cpymad_logfile)

