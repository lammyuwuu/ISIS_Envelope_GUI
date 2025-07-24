from helper_functions import *
from cpymad_helpers import *
from cpymad_closed_orbit_matching_functions import *
from ISIS_tune_control_functions import *




######################################################### 
#*******************  FUNCTIONS AND CLASSES FROM NOTEBOOK ***********************
#########################################################

def cpymad_set_isis_cycle_time(madx_instance, max_E, time):
    # Ensure time is a float and in valid increments
    if not isinstance(time, float) or time < 0.0 or time > 10.0 or (time * 10) % 5 != 0:
        print(f"Error: time must be a float between 0.0 and 10.0 in 0.5 increments. Received: {time}")
        return

    # Generate dataframe of synchrotron energy and related info
    energy_df = synchrotron_energy_df(max_E, intervals=20)

    # store some values for this time point
    try:
        energy = energy_df[energy_df['Time [ms]'] == time]['Energy [MeV]'].iloc[0]
        pc = energy_df[energy_df['Time [ms]'] == time]['Momentum [GeV/c]'].iloc[0]
    except IndexError:
        print(f"Error: No matching time value found in energy dataframe for time = {time} ms")
        return

    # set the beam to this energy in cpymad
    madx_instance.input(f'beam, particle = proton, pc = {pc};')

    # print confirmation
    print(f'ISIS cpymad run, energy set to {energy} MeV, pc = {pc}')


def set_DW_tunes(madx_instance, cpymad_logfile, input_df, sequence_name='synchrotron'):
    madx_instance.globals.kqtd = 0
    madx_instance.globals.kqtf = 0

    Iqtf, Iqtd = tune_to_trim_quad_current_di(Qh=input_df['Qh'], Qv=input_df['Qv'], baseQh=4.331, baseQv=3.731, pn=input_df['pn'], z=np.array([-4.73e-3, -5.99E-03, 4.45E-03, 2.40E-03]))
    Kqtf = current_to_strength(Iqtf, Gcal=1.997e-3, Brho=input_df['Rigidity'], pn=input_df['pn'])
    Kqtd = current_to_strength(Iqtd, Gcal=1.997e-3, Brho=input_df['Rigidity'], pn=input_df['pn'])
    madx_instance.globals.kqtd = Kqtd
    madx_instance.globals.kqtf = Kqtf
    
    twiss_df = cpymad_madx_twiss(madx_instance, cpymad_logfile, sequence_name)  
    return madx_instance.table.summ.q1[0], madx_instance.table.summ.q2[0]



def tune_calculation_iterator(madx_instance, df_input, cpymad_logfile, sequence_name='synchrotron'):
    qx_array = []  # to store output
    qy_array = []  # to store output
    set_qy_array = []
    qx_out = 0
    qy_out = 0 
    
    # Add new columns to df_input: 'Machine Qh', 'Machine Qv'
    df_input['Machine Qh'] = None
    df_input['Machine Qv'] = None
    
    for index, row in df_input.iterrows():
        
        # Generate new df of single row
        df_cut = row
        try: qx_out, qy_out = set_DW_tunes(madx_instance, cpymad_logfile, row, sequence_name)
        except: print('skipped point ', row['Qh'], row['Qv'])
        
        # Append value for Machine Qh = qx_out in df_input
        df_input.at[index, 'Machine Qh'] = qx_out
        # Append value for Machine Qv = qy_out in df_input
        df_input.at[index, 'Machine Qv'] = qy_out

    return df_input

class resonance_lines(object):

    def __init__(self, Qx_range, Qy_range, orders, periodicity, legend=False):

        if np.std(Qx_range):
            self.Qx_min = np.min(Qx_range)
            self.Qx_max = np.max(Qx_range)
        else:
            self.Qx_min = np.floor(Qx_range)-0.05
            self.Qx_max = np.floor(Qx_range)+1.05
        if np.std(Qy_range):
            self.Qy_min = np.min(Qy_range)
            self.Qy_max = np.max(Qy_range)
        else:
            self.Qy_min = np.floor(Qy_range)-0.05
            self.Qy_max = np.floor(Qy_range)+1.05

        self.periodicity = periodicity
        self.legend_flag = legend

        nx, ny = [], []

        for order in np.nditer(np.array(orders)):
            t = np.array(range(-order, order+1))
            nx.extend(order - np.abs(t))
            ny.extend(t)
        nx = np.array(nx)
        ny = np.array(ny)

        cextr = np.array([nx*np.floor(self.Qx_min)+ny*np.floor(self.Qy_min), \
                          nx*np.ceil(self.Qx_max)+ny*np.floor(self.Qy_min), \
                          nx*np.floor(self.Qx_min)+ny*np.ceil(self.Qy_max), \
                          nx*np.ceil(self.Qx_max)+ny*np.ceil(self.Qy_max)], dtype='int')
        cmin = np.min(cextr, axis=0)
        cmax = np.max(cextr, axis=0)
        res_sum = [range(cmin[i], cmax[i]+1) for i in range(cextr.shape[1])]
        self.resonance_list = zip(nx, ny, res_sum)

    def plot_resonance_ax(self, ax1):
        # Remove Borders
        #ax1.spines['top'].set_visible(False);
        #ax1.spines['bottom'].set_visible(False);
        #ax1.spines['left'].set_visible(False);
        #ax1.spines['right'].set_visible(False);
        #ax1.axes.get_yaxis().set_visible(False); 
        #ax1.axes.get_xaxis().set_visible(False);  
        
        Qx_min = self.Qx_min
        Qx_max = self.Qx_max
        Qy_min = self.Qy_min
        Qy_max = self.Qy_max 
        ax1.set_xlim(Qx_min, Qx_max);
        ax1.set_ylim(Qy_min, Qy_max);        
        ax1.set_xlabel(r'Horizontal Tune $Q_x$')
        ax1.set_ylabel(r'Vertical Tune $Q_y$')

        for resonance in self.resonance_list:
            nx = resonance[0]
            ny = resonance[1]
            for res_sum in resonance[2]:
                if ny:                    
                    
                    if ny%2:                                                  
                        if res_sum%self.periodicity:
                            ax1.plot([Qx_min, Qx_max], [(res_sum-nx*Qx_min)/ny, (res_sum-nx*Qx_max)/ny], ls='--', color='b', lw=0.5, label='Non-Systematic Skew')
                        else:
                            ax1.plot([Qx_min, Qx_max], [(res_sum-nx*Qx_min)/ny, (res_sum-nx*Qx_max)/ny], ls='--', color='r', lw=1, label='Systematic Skew')               
                        
                    else:
                        if res_sum%self.periodicity:
                            ax1.plot([Qx_min, Qx_max], [(res_sum-nx*Qx_min)/ny, (res_sum-nx*Qx_max)/ny], color='b', lw=0.5, label='Non-Systematic Normal')
                        else:
                            ax1.plot([Qx_min, Qx_max], [(res_sum-nx*Qx_min)/ny, (res_sum-nx*Qx_max)/ny], color='r', lw=1, label='Systematic Normal')
                else:                
                    if res_sum%self.periodicity:
                        ax1.plot([float(res_sum)/nx, float(res_sum)/nx],[Qy_min, Qy_max], color='b', lw=0.5, label='Non-Systematic Normal')
                    else:
                        ax1.plot([float(res_sum)/nx, float(res_sum)/nx],[Qy_min, Qy_max], color='r', lw=1, label='Systematic Normal')
    
        if self.legend_flag: 
            custom_lines = [Line2D([0], [0], color='r', lw=4),
                Line2D([0], [0], color='b', lw=4),
                Line2D([0], [2], color='r', lw=4, ls='--'),
                Line2D([0], [2], color='b', lw=4, ls='--')]
            ax1.legend(custom_lines, ['Systematic Normal', 'Non-Systematic Normal', 'Systematic Skew', 'Non-Systematic Skew'])
    
    
    def plot_resonance_fig(self, figure_object = None):
        plt.ion()
        if figure_object:
            fig = figure_object
            plt.figure(fig.number)
        else:
            fig = plt.figure()
        Qx_min = self.Qx_min
        Qx_max = self.Qx_max
        Qy_min = self.Qy_min
        Qy_max = self.Qy_max 
        plt.xlim(Qx_min, Qx_max)
        plt.ylim(Qy_min, Qy_max)
        plt.xlabel(r'Horizontal Tune $Q_x$')
        plt.ylabel(r'Vertical Tune $Q_y$')
        for resonance in self.resonance_list:
            nx = resonance[0]
            ny = resonance[1]
            for res_sum in resonance[2]:
                if ny:
                    line, = plt.plot([Qx_min, Qx_max], \
                        [(res_sum-nx*Qx_min)/ny, (res_sum-nx*Qx_max)/ny])
                else:
                    line, = plt.plot([float(res_sum)/nx, float(res_sum)/nx],[Qy_min, Qy_max])
                    
                    
                    
                if ny%2:
                    plt.setp(line, linestyle='--') # for skew resonances
                if res_sum%self.periodicity:
                    plt.setp(line, color='b') # non-systematic resonances
                else:
                    plt.setp(line, color='r', linewidth=2.0) # systematic resonances
        plt.draw()
        return fig
    
    def print_resonances(self):
        for resonance in self.resonance_list:
            for res_sum in resonance[2]:
                print_string =  '%s %s%s = %s\t%s'%(str(resonance[0]).rjust(2), ("+", "-")[resonance[1]<0], \
                        str(abs(resonance[1])).rjust(2), str(res_sum).rjust(4), \
                        ("(non-systematic)", "(systematic)")[res_sum%self.periodicity==0])
                print(print_string)

######################################################### 
#******************* END OF FUNCTIONS AND CLASSES FROM NOTEBOOK ***********************
#########################################################

def get_harmonic(pv_name: str, cycletime: str | float) -> float:
    """Return rows matching a specific CycleTime and PV value."""
    # Get the directory of the current script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    csv_path = os.path.join(script_dir, "..", "Collected_EPICS_data", "get_EPICS_Harmonics_full_cycle.dat")

    harmonic_data = pd.read_csv(csv_path)

    harmonic_index = harmonic_data["PV"].str.contains(f"DWTRIM::{pv_name}:AT_TIME:{cycletime}MS")
    return harmonic_data[harmonic_index]["Harmonic"].to_list()[0]



def get_twiss_table(time_point, has_harmonic):
    # 1. Instantiate a MAD-X object
    cpymad_logfile = 'cpymad_logfile.txt'
    sequence_name = 'synchrotron'

    madx = cpymad_start(cpymad_logfile)

    #`madx` is the simulation object, next we have to give it the description of the synchrotron, also called the **lattice**
    #We have a number of lattices to choose from

    # lattice_folder = '../Lattice_Files/00_Simplified_Lattice/'
    # lattice_folder = '../Lattice_Files/01_Original_Lattice/'
    # lattice_folder = '../Lattice_Files/02_Aperture_Lattice/'
    # lattice_folder = '../Lattice_Files/03_CO_Kick_Lattice/'
    lattice_folder = '../Lattice_Files/04_New_Harmonics/'

    madx.call(file=lattice_folder+'ISIS.injected_beam')
    madx.call(file=lattice_folder+'ISIS.strength')
    madx.call(file=lattice_folder+'2023.strength')
    madx.call(file=lattice_folder+'ISIS.elements')
    madx.call(file=lattice_folder+'ISIS.sequence')

    cpymad_check_and_use_sequence(madx, cpymad_logfile, sequence_name)

    ## Set harmonics

    if has_harmonic:
        madx.globals['D7COS'] = get_harmonic("D7COS", time_point)
        madx.globals['D8COS'] = get_harmonic("D8COS", time_point)
        madx.globals['F8COS'] = get_harmonic("F8COS", time_point)
        madx.globals['D7SIN'] = get_harmonic("D7SIN", time_point)
        madx.globals['D8SIN'] = get_harmonic("D8SIN", time_point)
        madx.globals['F8SIN'] = get_harmonic("F8SIN", time_point)



    # 2. Set the cycle time = beam energy
    max_E = 800 # 800 MeV
    cycle_time = time_point # 1.5 milliseconds into the 10 ms acceleration cycle of the synchrotron
    cpymad_set_isis_cycle_time(madx, max_E, cycle_time)

    twiss_0 = cpymad_madx_twiss(madx, cpymad_logfile, sequence_name)




    # 3. Calculate lattice parameters using **TWISS** command



    cpymad_plot_madx_twiss_quads(madx, twiss_0, title='Initial lattice tune')
    return twiss_0


