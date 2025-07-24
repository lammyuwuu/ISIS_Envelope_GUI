#!/usr/bin/env python
# cpymad functions to be used with modified ISIS lattice for closed orbit calculations
# Created: 17.06.22 Haroon Rafique STFC ISIS Accelerator Division
# Requires cpymad tfs-pandas python libraries

from helper_functions import *
from cpymad_helpers import *

def set_tune_DW(madx_instance, cpymad_logfile, Qh, Qv, time, sequence_name='synchrotron'):
    tq_currents_df = tune_di_df(Qh=np.array([Qh]), Qv=np.array([Qv]), baseQh=4.331, baseQv=3.731, time_array=np.array([time]), z=np.array([-4.73e-3, -5.99E-03, 4.45E-03, 2.40E-03]))
    
    madx_instance.globals.kqtd = tq_currents_df['K_qtd'][0]
    madx_instance.globals.kqtf = tq_currents_df['K_qtf'][0]
    
    twiss_df = cpymad_madx_twiss(madx_instance, cpymad_logfile, sequence_name)  
    return madx_instance.table.summ.q1[0], madx_instance.table.summ.q2[0]
    
def get_isis_tunes(madx_instance, cpymad_logfile):
    initial_tunes = {}
    
    # Helper function for repeated tasks
    def get_tune_data(madx, cpymad_log, seq_name, is_ptc=False):
        if is_ptc:
            cpymad_ptc_twiss(madx, cpymad_log, seq_name)
            qx, qy = madx.table.ptc_twiss_summary.q1[0], madx.table.ptc_twiss_summary.q2[0]
        else:
            cpymad_madx_twiss(madx, cpymad_log, seq_name)
            qx, qy = madx.table.summ.q1[0], madx.table.summ.q2[0]
            
        return qx, qy
    
    # Fetch tunes from superperiod matrix phase advance
    qx, qy = superperiod_matrix_phase_advance(madx_instance, cpymad_logfile, madx_instance.sequence.superperiod)
    initial_tunes['SP Matrix Qx'], initial_tunes['SP Matrix Qy'] = round(qx * 10, 4), round(qy * 10, 4)
    
    # Define list of operation settings for looping
    operations = [
        ('superperiod', 'MAD-X sp Qx', 'MAD-X sp Qy', False),
        ('superperiod', 'PTC sp Qx', 'PTC sp Qy', True),
        ('synchrotron', 'MAD-X Qx', 'MAD-X Qy', False),
        ('synchrotron', 'PTC Qx', 'PTC Qy', True)
    ]

    # Get tune data for each operation
    for seq_name, qx_key, qy_key, is_ptc in operations:
        qx, qy = get_tune_data(madx_instance, cpymad_logfile, seq_name, is_ptc)
        
        if 'sp' in qx_key:
            qx *= 10
            qy *= 10
        elif qx_key == 'PTC Qx':
            qx += 4
            qy += 3
        
        # Rounding the values to 4 decimal places
        initial_tunes[qx_key], initial_tunes[qy_key] = round(qx, 4), round(qy, 4)
    
    return initial_tunes
        
def isis_reset_trim_quads(madx_instance):   
    
    # trim quads
    madx_instance.globals.kqtd_0 = 'kqtd + HER0qtd'
    madx_instance.globals.kqtd_1 = 'kqtd + HER1qtd'
    madx_instance.globals.kqtd_2 = 'kqtd + HER2qtd'
    madx_instance.globals.kqtd_3 = 'kqtd + HER3qtd'
    madx_instance.globals.kqtd_4 = 'kqtd + HER4qtd'
    madx_instance.globals.kqtd_5 = 'kqtd + HER5qtd'
    madx_instance.globals.kqtd_6 = 'kqtd + HER6qtd'     
    madx_instance.globals.kqtd_7 = 'kqtd + HER7qtd'
    madx_instance.globals.kqtd_8 = 'kqtd + HER8qtd'
    madx_instance.globals.kqtd_9 = 'kqtd + HER9qtd'
    
    madx_instance.globals.kqtf_0 = 'kqtf + HER0qtf'
    madx_instance.globals.kqtf_1 = 'kqtf + HER1qtf'
    madx_instance.globals.kqtf_2 = 'kqtf + HER2qtf'
    madx_instance.globals.kqtf_3 = 'kqtf + HER3qtf'
    madx_instance.globals.kqtf_4 = 'kqtf + HER4qtf'
    madx_instance.globals.kqtf_5 = 'kqtf + HER5qtf'
    madx_instance.globals.kqtf_6 = 'kqtf + HER6qtf'     
    madx_instance.globals.kqtf_7 = 'kqtf + HER7qtf'
    madx_instance.globals.kqtf_8 = 'kqtf + HER8qtf'
    madx_instance.globals.kqtf_9 = 'kqtf + HER9qtf'
    
def isis_print_trim_quads(madx_instance):   
    
    # Steering magnets
    print(madx_instance.globals.kqtd_0,
    madx_instance.globals.kqtd_1,
    madx_instance.globals.kqtd_2,
    madx_instance.globals.kqtd_3,
    madx_instance.globals.kqtd_4,
    madx_instance.globals.kqtd_5,
    madx_instance.globals.kqtd_6,      
    madx_instance.globals.kqtd_7,
    madx_instance.globals.kqtd_8,
    madx_instance.globals.kqtd_9,    
    madx_instance.globals.kqtf_0,
    madx_instance.globals.kqtf_1,
    madx_instance.globals.kqtf_2,
    madx_instance.globals.kqtf_3,
    madx_instance.globals.kqtf_4,
    madx_instance.globals.kqtf_5,
    madx_instance.globals.kqtf_6,      
    madx_instance.globals.kqtf_7,
    madx_instance.globals.kqtf_8,
    madx_instance.globals.kqtf_9)
    
def tune_to_trim_quad_current_di(Qh=4.331, Qv=3.731,
                                 baseQh=4.331, baseQv=3.731, pn=1.0,
                                 z=np.array([-4.73e-3, -5.99E-03, 4.45E-03, 2.40E-03])):
    '''
    Calculates the trim quad currents required to obtain the given tune under
    Di's tune control system. The default values for the base tunes and
    normalised momentum are those at t = 0.0 ms.

    Parameters
    ----------
    Qh : Float, optional
        Horizontal tune required.
        Can be a float array of tunes. Must be of same length as Qv.
        The default is 4.331.
    Qv : Float, optional
        Vertical tune required.
        Can be a float array of tunes. Must be of same length as Qh
        The default is 3.731.
    baseQh : Float, optional
        Base horizontal tune.
        Can be a float array of tunes. Must be of same length as Qh, Qv.
        The default is 4.331.
    baseQv : Float, optional
        Base vertical tune.
        Can be a float array of tunes. Must be of same length as Qh, Qv.
        The default is 3.731.
    pn : Float, optional
        Normalised momentum.
        Can be a float array of momenta. Must be of same length as Qh, Qv.
        The default is 1.0.
    z : Float array of length 4, optional
        Coefficients of the tune control system.
        The default is np.array([-4.73e-3, -5.99E-03, 4.45E-03, 2.40E-03]).

    Returns
    -------
    I_F : Float, or array of floats with length same as Qh, Qv
        Current to apply to QTF to obtain Qh, Qv according to Di controls.
    I_D : Float, or array of floats with length same as Qh, Qv
        Current to apply to QTD to obtain Qh, Qv according to Di controls.

    '''

    # Get the control coefficients
    z1, z2, z3, z4 = z

    # Calculate the change in tune required
    dQh = Qh - baseQh
    dQv = Qv - baseQv

    # Calculate the currents required
    I_F = pn * (z1 * dQv - z3 * dQh) / (z1 * z4 - z2 * z3)
    I_D = pn * (z4 * dQh - z2 * dQv) / (z1 * z4 - z2 * z3)

    # Return the currents
    return -I_F, -I_D    
    
def current_to_strength(I, Gcal=1.997e-3, Brho=1.23, pn=1.0):
    k = I * Gcal / Brho / pn
    return k
    
def tune_di_df_srm(Qh=4.331, Qv=3.731, baseQh=4.331, baseQv=3.731, time_array=None, E_Max=800, z=np.array([-4.73e-3, -5.99E-03, 4.45E-03, 2.40E-03])):
    if time_array is None: 
        print('tune_to_trim_quad_current_di_df::error: time_array is None')
        return False    
    else:
        pn_array = return_pn(time_array)
        df_KE = synchrotron_kinetic_energy_df(E_Max, time_array)
        
    Iqtf_array = []
    Iqtd_array = []
    Kqtf_array = []
    Kqtd_array = []
    
    for pn, (_, ke_row) in zip(pn_array, df_KE.iterrows()):
        Iqtf, Iqtd = tune_to_trim_quad_current_di(Qh, Qv, baseQh, baseQv, pn, z)
        Iqtf_array.append(Iqtf)
        Iqtd_array.append(Iqtd)
        Kqtf_array.append(current_to_strength(Iqtf, Gcal=1.997e-3, Brho=ke_row['Rigidity [Tm]'], pn=pn))
        Kqtd_array.append(current_to_strength(Iqtd, Gcal=1.997e-3, Brho=ke_row['Rigidity [Tm]'], pn=pn))
        
    data = {
        'time': time_array,
        'pn': pn_array,
        'I_qtf': Iqtf_array,
        'I_qtd': Iqtd_array,
        'K_qtf': Kqtf_array,
        'K_qtd': Kqtd_array
    }
    
    df = pd.DataFrame(data)
    return df
    
def tune_di_df(Qh, Qv, baseQh=4.331, baseQv=3.731, time_array=None, E_Max=800, z=np.array([-4.73e-3, -5.99E-03, 4.45E-03, 2.40E-03])):
    if time_array is None: 
        print('tune_to_trim_quad_current_di_df::error: time_array is None')
        return False    
    if len(Qh) != len(time_array) or len(Qv) != len(time_array):
        print('tune_to_trim_quad_current_di_df::error: Qh and Qv must have the same length as time_array')
        return False
    
    pn_array = return_pn(time_array)
    df_KE = synchrotron_kinetic_energy_df(E_Max, time_array)
        
    Iqtf_array = []
    Iqtd_array = []
    Kqtf_array = []
    Kqtd_array = []
    brho_array = []
    
    for qh, qv, pn, (_, ke_row) in zip(Qh, Qv, pn_array, df_KE.iterrows()):
        Iqtf, Iqtd = tune_to_trim_quad_current_di(qh, qv, baseQh, baseQv, pn, z)
        Iqtf_array.append(Iqtf)
        Iqtd_array.append(Iqtd)
        Kqtf_array.append(current_to_strength(Iqtf, Gcal=1.997e-3, Brho=ke_row['Rigidity [Tm]'], pn=pn))
        Kqtd_array.append(current_to_strength(Iqtd, Gcal=1.997e-3, Brho=ke_row['Rigidity [Tm]'], pn=pn))
        brho_array.append(ke_row['Rigidity [Tm]'])
        
    data = {
        'time': time_array,
        'pn': pn_array,
        'Qh': Qh,
        'Qv': Qv,
        'I_qtf': Iqtf_array,
        'I_qtd': Iqtd_array,
        'K_qtf': Kqtf_array,
        'K_qtd': Kqtd_array,
        'Rigidity': brho_array
    }
    
    df = pd.DataFrame(data)
    return df

    
def tune_to_phase_advance(Qy, superperiods=10, arccos_sign=1, arccos_period=0):
    phiy = 2 * np.cos(-arccos_sign * 2 * np.pi * (Qy / superperiods - arccos_period))
    return phiy
    
def sp_trace_coefficients(A, B, L_trim=0.307):
#def sp_trace_coefficients(B, A, L_trim=0.307): #PGH
    
    #c1 = round_sig(L_trim**2 * A[0,1] * B[0,1],5)
    #c2 = round_sig(L_trim*(A[0,0]*B[0,1]+A[0,1]*B[1,1]),5)
    #c3 = round_sig(L_trim*(A[0,1]*B[0,0]+A[1,1]*B[0,1]),5)
    #c4 = round_sig(A[0,0]*B[0,0] + A[0,1]*B[1,0] + A[1,0]*B[0,1] + A[1,1]*B[1,1],5)
    
    c1 = (L_trim**2 * A[0,1] * B[0,1])
    c2 = (L_trim*(A[0,0]*B[0,1]+A[0,1]*B[1,1]))
    c3 = (L_trim*(A[0,1]*B[0,0]+A[1,1]*B[0,1]))
    c4 = (A[0,0]*B[0,0] + A[0,1]*B[1,0] + A[1,0]*B[0,1] + A[1,1]*B[1,1])
    
    return [c1,c2,c3,c4]   
    
def phi_to_k_coeffs(phih, phiv, h, v):
    # Get one turn transfer matrix trace coefficients
    h1, h2, h3, h4 = h
    v1, v2, v3, v4 = v

    # Calculate the coefficients of the quadratic equation that gives k_QTF
    a = (-h1 * v3 - h3 * v1)
    b = (phih * v1 - phiv * h1) + ((h1 * v4 - h4 * v1) + (h3 * v2 - h2 * v3))
    c = (-phih * v2 - phiv * h2) + (h2 * v4 + h4 * v2)

    # Return the coefficients
    return a, b, c
    
def find_min_k(k_F1, k_F2, k_D1, k_D2):
    k_F, k_D = np.where((abs(k_F1) + abs(k_D1)) < (abs(k_F2) + abs(k_D2)),
                        np.array([k_F1, k_D1]), np.array([k_F2, k_D2]))   
    return k_F, k_D
    
def coeffs_to_strengths(a, b, c, phih, h):
    # Get one turn transfer matrix trace coefficients
    h1, h2, h3, h4 = h

    # Solve the simultaneous equations for k_QTF, k_QTD
    k_D1 = (-b + (b**2 - 4 * a * c)**0.5) / (2 * a)
    k_D2 = (-b - (b**2 - 4 * a * c)**0.5) / (2 * a)
    k_F1 = (phih - h3 * k_D1 - h4) / (h1 * k_D1 + h2)
    k_F2 = (phih - h3 * k_D2 - h4) / (h1 * k_D2 + h2)

    # Find the solutions with the minimum trim quad strengths
    k_F, k_D = find_min_k(k_F1, k_F2, k_D1, k_D2)

    # Return trim quad strengths
    return k_F, k_D
    
def stable_sp_check(phix, phiy):
    
    x_stable = np.where(np.abs(phix)<=2, True, False)
    y_stable = np.where(np.abs(phiy)<=2, True, False)

    return np.logical_and(x_stable, y_stable)    
    
def test_sp(sp):       
    tolerance = 1E-4
    if np.abs(np.linalg.det(sp) - 1) >= tolerance:
        return False
    else:
        return matrix_phase_advance(sp)
        
def thin_quad_matrix_6(K,L):
    out_matrix = np.matrix(np.zeros(shape=[6,6]))    
    m_k = thin_quad_matrix_2(K, L)
    
    out_matrix[0,0] = m_k[0,0]
    out_matrix[0,1] = m_k[0,1]
    out_matrix[1,0] = m_k[1,0]
    out_matrix[1,1] = m_k[1,1]
    
    out_matrix[2,2] = m_k[0,0]
    out_matrix[2,3] = m_k[0,1]
    out_matrix[3,2] = -m_k[1,0]
    out_matrix[3,3] = m_k[1,1]
    
    out_matrix[4,4] = 1.0
    out_matrix[5,5] = 1.0
    
    return out_matrix
    
def thin_quad_matrix_2(K,L):
    return np.array([[1,0],[-K*L, 1]])
        
# Note that the modified lattice files for tune controls must be used - the one that has sequences sp_tune_m0 and sp_tune_m1
def get_lattice_tune_control_parameters(lat = './Lattice/', save_folder = 'cpymad_tune_control', scaling = None, requested_qx=4.31, requested_qy=3.79, verbose=False):
        
    if save_folder[-1] != '/': save_folder = save_folder+'/'
  
    sp_checks = False
    
    m0_seq = 'sp_tune_m0'
    m1_seq = 'sp_tune_m1'
    
    # start cpymad run
    logfile = 'cpymad_logfile_'+str(requested_qx)+'_'+str(requested_qy)+'.log'
    (madx_0, cpymad_logfile) = cpymad_start_ISIS(lat, save_folder, m0_seq, False, logfile)
    
    # Choose scaled lattice strengths (scaled quad K1 in main dipole, main quads, fringes, etc to match bare tune according to PGH)
    if scaling == None:
        # no scaling = (0.4336, 0.383)
        if verbose: print('get_lattice_tune_control_parameters: No lattice scaling applied')
        
    elif scaling == 'q':
        madx_0.call(file=lat+'Scaled_quad_fringes.strength') # (0.4323, 0.3763)
        if verbose: print('get_lattice_tune_control_parameters: scaling applied including quad fringes')
        
    elif scaling == 'n':
        madx_0.call(file=lat+'Scaled_no_fringes.strength') # (0.4312, 0.3789)
        if verbose: print('get_lattice_tune_control_parameters: scaling applied not including fringes')
        
    elif scaling == ('qd' or 'dq'):
        madx_0.call(file=lat+'Scaled_both_fringes.strength') # (0.4312, 0.3788)
        if verbose: print('get_lattice_tune_control_parameters: scaling applied including dipole and quad fringes')
    
    else:
        raise ValueError('get_lattice_tune_control_parameters: scaling not recognised, please select:\n\'q\' for scaling including quadrupole fringes \n\'n\' for scaling with no fringes\n\'qd\' for scaling with both quad and dipole fringes')
        
    # First Twiss    
    twiss_sp_0_file = str(save_folder+'twiss_sp_tune_m0.tfs')
    twiss_sp_0 = cpymad_madx_twiss(madx_0, cpymad_logfile, m0_seq, twiss_sp_0_file)
    # cpymad_plot_madx_twiss(madx_0, twiss_sp_0, title='ISIS sp_tune_m0', savename=(save_folder +'Twiss_sp_tune_m0.png'))  
        
    # El list 1
    sp_elements_0, superperiod_seq  = cpymad_sequence_to_dict(madx_0, madx_0.sequence.sp_tune_m0)
        
    if sp_checks:
        # Check tune 1
        sp0 = transport_matrix_between_elements(madx_0, sp_elements_0, 'sp_tune_m0$start', 'sp_tune_m0$end')    
        check1 = test_sp(sp0)    
        if check1 and verbose: print (check1)
        else: raise ValueError('get_lattice_tune_control_parameters: sp_tune_m0 determinant check failed')
            
    # Matrix M0
    m_0 = transport_matrix_between_elements(madx_0, sp_elements_0, 'sp_tune_m0$start', 'sp_dqtd2')
        
    # Next sequence and Twiss
    twiss_sp_1_file = str(save_folder+'twiss_sp_tune_m1.tfs')
    twiss_sp_1 = cpymad_madx_twiss(madx_0, cpymad_logfile, m1_seq, twiss_sp_1_file)
    # cpymad_plot_madx_twiss(madx_0, twiss_sp_1, title='ISIS sp_tune_m1', savename=(save_folder +'Twiss_sp_tune_m1.png'))  

    # El list 2
    sp_elements_1, superperiod_seq  = cpymad_sequence_to_dict(madx_0, madx_0.sequence.sp_tune_m1)
    
    if sp_checks:   
        # Check tune 2
        sp1 = transport_matrix_between_elements(madx_0, sp_elements_1, 'sp_tune_m1$start', 'sp_tune_m1$end') 
        check2 = test_sp(sp1)    
        if check2 and verbose: print (check2)
        else: raise ValueError('get_lattice_tune_control_parameters: sp_tune_m1 determinant check failed')
            
    # Matrix M1
    m_1 = transport_matrix_between_elements(madx_0, sp_elements_1, 'sp_tune_m1$start', 'sp_dqtf')

    if sp_checks:       
        # Check matrix sum
        m_2 = m_0 @ m_1
        check3 = test_sp(m_2)    
        if check3 and verbose: print (check3)
        else: raise ValueError('get_lattice_tune_control_parameters: sum matrix m_2 determinant check failed')        

    # decouple transverse matrices
    m_0_h, m_0_v = helper_functions.transport_matrix_transverse_decomposition(m_0)
    m_1_h, m_1_v = helper_functions.transport_matrix_transverse_decomposition(m_1)

    # h and v coefficients
    # Note HR M_1 == PGH A, HR M_0 == PGH B matrices
    h_c = sp_trace_coefficients(m_1_h, m_0_h)
    v_c = sp_trace_coefficients(m_1_v, m_0_v)
    
    # Convert requested tunes to phase advance
    phi_h = tune_to_phase_advance(requested_qx)
    phi_v = tune_to_phase_advance(requested_qy)
    
    # Stability check
    if not stable_sp_check(phi_h, phi_v): raise ValueError('get_lattice_tune_control_parameters: stable_sp_check failed: requested tunes outside area of stability')
    
    # Extract quadratic coefficients
    a, b, c = phi_to_k_coeffs(phi_h, phi_v, h_c, v_c)
    
    # Convert to quadrupole strengths for MAD input
    k_F, k_D = coeffs_to_strengths(a, b, c, phi_h, h_c)
        
    if sp_checks and verbose: 
        print('M0 Tune = ', check1)
        print('M1 Tune = ', check2)
        print('M2 Tune = ', check3)
        
    # close cpymad run
    madx_0.quit()
        
    # Difference between PGH Tracker and cpymad ???
    # Result: K signs are opposite
    #return k_F, k_D
    return -k_F, -k_D

def calculate_tune_error(lat = './Lattice/', save_folder = 'cpymad_tune_control/', scaling = None, requested_qx=4.31, requested_qy=3.79, k_F = 0.0, k_D = 0.0, verbose=False):

    if save_folder[-1] != '/': save_folder = save_folder+'/'
    
    if k_D == k_F == 0.0: raise ValueError('calculate_tune_error: k_D and k_F set to zero')
    
    sequence='synchrotron'
    
    # start cpymad run
    (madx_0, cpymad_logfile) = cpymad_start_ISIS(lat, save_folder, sequence)
    
    # Choose scaled lattice strengths (scaled quad K1 in main dipole, main quads, fringes, etc to match bare tune according to PGH)
    if scaling == None:
        # no scaling = (0.4336, 0.383)
        if verbose: print('get_lattice_tune_control_parameters: No lattice scaling applied')
        
    elif scaling == 'q':
        madx_0.call(file=lat+'Scaled_quad_fringes.strength') # (0.4323, 0.3763)
        if verbose: print('calculate_tune_error: scaling applied including quad fringes')
        
    elif scaling == 'n':
        madx_0.call(file=lat+'Scaled_no_fringes.strength') # (0.4312, 0.3789)
        if verbose: print('calculate_tune_error: scaling applied not including fringes')
        
    elif scaling == ('qd' or 'dq'):
        madx_0.call(file=lat+'Scaled_both_fringes.strength') # (0.4312, 0.3788)
        if verbose: print('calculate_tune_error: scaling applied including dipole and quad fringes')
    
    else:
        raise ValueError('calculate_tune_error: scaling not recognised, please select:\n\'q\' for scaling including quadrupole fringes \n\'n\' for scaling with no fringes\n\'qd\' for scaling with both quad and dipole fringes')
    
    # Set the trim qaud strengths
    madx_0.globals.kqtf = k_F
    madx_0.globals.kqtd = k_D        
    
    # Twiss to get tunes
    twiss_0_file = str(save_folder+'twiss_tune_test_'+str()+'_'+str()+'.tfs')
    twiss_0 = cpymad_madx_twiss(madx_0, cpymad_logfile, sequence, twiss_0_file)
    twiss_plot_title = str('ISIS Tune Test: Requested (Qx, Qy)='+str(requested_qx)+','+str(requested_qy))
    twiss_plot_file = str(save_folder+'twiss_tune_test_'+str(requested_qx)+'_'+str(requested_qy)+'.png')
    cpymad_plot_madx_twiss(madx_0, twiss_0, title=twiss_plot_title, savename=twiss_plot_file)  
        
    Qx = madx_0.table.summ.q1[0]
    Qy = madx_0.table.summ.q2[0]
    
    abs_qx = np.abs(requested_qx - Qx)
    abs_qy = np.abs(requested_qy - Qy)
        
    Q = [Qx, Qy]
    Q_err = [abs_qx, abs_qy]
        
    # close cpymad run
    madx_0.quit()
    
    return Q, Q_err
    
def calculate_tune_error_sm(lat = './Lattice/', save_folder = 'cpymad_tune_control/', scaling = None, requested_qx=4.31, requested_qy=3.79, k_F = 0.0, k_D = 0.0, verbose=False):
    
    if save_folder[-1] != '/': save_folder = save_folder+'/'
    
    if k_D == k_F == 0.0: raise ValueError('calculate_tune_error_sm: k_D and k_F set to zero')
    
    m0_seq = 'sp_tune_m0'
    m1_seq = 'sp_tune_m1'
    
    # start cpymad run
    (madx_0, cpymad_logfile) = cpymad_start_ISIS(lat, save_folder, m0_seq)
    
    # Choose scaled lattice strengths (scaled quad K1 in main dipole, main quads, fringes, etc to match bare tune according to PGH)
    if scaling == None:
        # no scaling = (0.4336, 0.383)
        if verbose: print('calculate_tune_error_sm: No lattice scaling applied')
        
    elif scaling == 'q':
        madx_0.call(file=lat+'Scaled_quad_fringes.strength') # (0.4323, 0.3763)
        if verbose: print('calculate_tune_error_sm: scaling applied including quad fringes')
        
    elif scaling == 'n':
        madx_0.call(file=lat+'Scaled_no_fringes.strength') # (0.4312, 0.3789)
        if verbose: print('calculate_tune_error_sm: scaling applied not including fringes')
        
    elif scaling == ('qd' or 'dq'):
        madx_0.call(file=lat+'Scaled_both_fringes.strength') # (0.4312, 0.3788)
        if verbose: print('calculate_tune_error_sm: scaling applied including dipole and quad fringes')
    else:
        raise ValueError('calculate_tune_error_sm: scaling not recognised, please select:\n\'q\' for scaling including quadrupole fringes \n\'n\' for scaling with no fringes\n\'qd\' for scaling with both quad and dipole fringes')
        
    # First Twiss    
    twiss_sp_0_file = str(save_folder+'twiss_sp_tune_m0.tfs')
    twiss_sp_0 = cpymad_madx_twiss(madx_0, cpymad_logfile, m0_seq, twiss_sp_0_file)
    #cpymad_plot_madx_twiss(madx_0, twiss_sp_0, title='ISIS sp_tune_m0', savename=(save_folder +'Twiss_sp_tune_m0.png'))  
        
    # El list 1
    sp_elements_0, superperiod_seq  = cpymad_sequence_to_dict(madx_0, madx_0.sequence.sp_tune_m0)
        
    # Check tune 1
    sp0 = transport_matrix_between_elements(madx_0, sp_elements_0, 'sp_tune_m0$start', 'sp_tune_m0$end')    
    #check1 = test_sp(sp0)    
    #if check1 and verbose: print (check1)
    #else: raise ValueError('calculate_tune_error_sm: sp_tune_m0 determinant check failed')
        
    # Matrix M0
    m_0 = transport_matrix_between_elements(madx_0, sp_elements_0, 'sp_dqtf2', 'sp_dqtd2')
        
    # Next sequence and Twiss
    twiss_sp_1_file = str(save_folder+'twiss_sp_tune_m1.tfs')
    twiss_sp_1 = cpymad_madx_twiss(madx_0, cpymad_logfile, m1_seq, twiss_sp_1_file)
    #cpymad_plot_madx_twiss(madx_0, twiss_sp_1, title='ISIS sp_tune_m1', savename=(save_folder +'Twiss_sp_tune_m1.png'))  

    # El list 2
    sp_elements_1, superperiod_seq  = cpymad_sequence_to_dict(madx_0, madx_0.sequence.sp_tune_m1)
    
    # Check tune 2
    sp1 = transport_matrix_between_elements(madx_0, sp_elements_1, 'sp_tune_m1$start', 'sp_tune_m1$end') 
    #check2 = test_sp(sp1)    
    #if check2 and verbose: print (check2)
    #else: raise ValueError('calculate_tune_error_sm: sp_tune_m1 determinant check failed')
        
    # Matrix M1
    m_1 = transport_matrix_between_elements(madx_0, sp_elements_1, 'sp_dqtd', 'sp_dqtf')
    
    # Check matrix sum
    m_2 = m_0 @ m_1
    #check3 = test_sp(m_2)    
    #if check3 and verbose: print (check3)
    #else: raise ValueError('calculate_tune_error_sm: sum matrix m_2 determinant check failed')

    # Calculate tune from SP matrix  
    
    # Need thin lens matrices for trim quads
    m_kf = thin_quad_matrix_6(k_F, 0.307)
    m_kd = thin_quad_matrix_6(k_D, 0.307)
    
    # multiply matrices
    matrix_dict = {}
    matrix_dict['qtd'] = m_kd
    matrix_dict['m1'] = m_1
    matrix_dict['qtf'] = m_kf
    matrix_dict['m0'] = m_0
    
    m_sp = multiply_matrices(matrix_dict)
    check4 = test_sp(m_sp)   
    if check4 and verbose: print (check4)
        
    # Get tunes
    qx, qy = matrix_phase_advance(m_sp)
    
    abs_qx = np.abs(requested_qx - qx*10)
    abs_qy = np.abs(requested_qy - qy*10)
        
    Q = [qx*10, qy*10]
    Q_err = [abs_qx, abs_qy]
        
    # close cpymad run
    madx_0.quit()
    
    return Q, Q_err
