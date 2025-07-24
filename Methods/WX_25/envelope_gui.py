import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import pandas as pd
from datetime import datetime

# Import necessary classes and functions
from ISIS_tune_control_functions import *
from cpymad_closed_orbit_matching_functions import *
from helper_functions import round_sig, synchrotron_energy_data, synchrotron_energy_df
from cpymad_helpers import cpymad_madx_twiss, cpymad_start, cpymad_check_and_use_sequence
from get_tune_values import *
from plotly.subplots import make_subplots
from plot_tune import *

cpymad_logfile = 'cpymad_logfile.txt'
sequence_name = 'synchrotron'
madx = cpymad_start(cpymad_logfile)
lattice_folder = '../../Lattice_Files/00_Simplified_Lattice/'
madx.call(file=lattice_folder+'ISIS.injected_beam')
madx.call(file=lattice_folder+'ISIS.strength')
madx.call(file=lattice_folder+'2023.strength')
madx.call(file=lattice_folder+'ISIS.elements')
madx.call(file=lattice_folder+'ISIS.sequence')
cpymad_check_and_use_sequence(madx, cpymad_logfile, sequence_name)

class cpymad_ErrorTableBuilder:
    def __init__(self, twiss_df):
        """
        Initialize the error table builder with a MAD-X twiss table DataFrame.
        """
        self.twiss_df = twiss_df.copy()
        self.error_df = pd.DataFrame(columns=["NAME", "DX", "DY", "DS", "DPHI", "DTHETA", "DPSI", "S"])

    def match_main_magnet_parts(self, base_name):
        """
        Return a list of element names matching base_name (e.g. 'sp6_qf'), 
        excluding fringe-like variants and stripping colon suffixes.
        """
        names = self.twiss_df['name'].astype(str)
        matched = names[
            names.str.startswith(base_name) &
            names.str.len().ge(len(base_name)) &
            names.str[len(base_name)].isin(['_', ':'])
        ]
        return matched.str.split(':').str[0].unique().tolist()      

    def add_dipole_misalignment(self, base_name, misalignment_type, value_mm):
        """
        Add a misalignment value (in mm or mrad) for all parts of a dipole element.
    
        Parameters:
            - base_name: base name of dipole (e.g. "sp6_dip")
            - misalignment_type: one of "DX", "DY", "DS", "DPHI", "DTHETA", "DPSI"
            - value_mm: value in mm or mrad (e.g. 1.5 for 1.5 mm)
        """
        misalignment_type = misalignment_type.upper()
        assert misalignment_type in ["DX", "DY", "DS", "DPHI", "DTHETA", "DPSI"], f"Invalid type {misalignment_type}"
    
        names = self.twiss_df['name'].astype(str)
        matched = names[names.str.startswith(base_name)]
        parts = matched.str.split(":").str[0].unique().tolist()
    
        # Prepare Twiss lookup with cleaned names (without colon suffixes)
        twiss_lookup = (
            self.twiss_df.copy()
            .assign(name_clean=self.twiss_df["name"].astype(str).str.split(":").str[0])
            .drop_duplicates(subset="name_clean")
            .set_index("name_clean")["s"]
        )
    
        for name in parts:
            match = self.error_df["NAME"] == name
            if match.any():
                self.error_df.loc[match, misalignment_type] = float(value_mm)
            else:
                try:
                    s_val = float(twiss_lookup.get(name))
                except Exception:
                    s_val = float("inf")
    
                columns = ["NAME", "DX", "DY", "DS", "DPHI", "DTHETA", "DPSI", "S"]
                row_full = {col: 0.0 for col in columns}
                row_full.update({
                    "NAME": name,
                    misalignment_type: float(value_mm),
                    "S": s_val
                })
                self.error_df = pd.concat([self.error_df, pd.DataFrame([row_full])], ignore_index=True)

        
    def add_quadrupole_misalignment(self, base_name, misalignment_type, value_mm):
        """
        Add a misalignment value (in mm or mrad) for all parts of a given magnet.
        
        Parameters:
            - base_name: base name of magnet (e.g. "sp6_qf")
            - misalignment_type: one of "DX", "DY", "DS", "DPHI", "DTHETA", "DPSI"
            - value_mm: value in mm or mrad (e.g. 1.5 for 1.5 mm)
        """
        misalignment_type = misalignment_type.upper()
        assert misalignment_type in ["DX", "DY", "DS", "DPHI", "DTHETA", "DPSI"], f"Invalid type {misalignment_type}"
        parts = self.match_main_magnet_parts(base_name)
        
        # Prepare Twiss lookup with cleaned names (without colon suffixes)
        twiss_lookup = (
            self.twiss_df.copy()
            .assign(name_clean=self.twiss_df["name"].astype(str).str.split(":").str[0])
            .drop_duplicates(subset="name_clean")
            .set_index("name_clean")["s"]
        )
    
        for name in parts:
            match = self.error_df["NAME"] == name
            if match.any():
                # Update existing row
                self.error_df.loc[match, misalignment_type] = float(value_mm)
            else:
                # Add new row
                try:
                    s_val = float(twiss_lookup.get(name))
                except Exception:
                    s_val = float("inf")
    
                row = {
                    "NAME": name,
                    "DX": 0.0,
                    "DY": 0.0,
                    "DS": 0.0,
                    "DPHI": 0.0,
                    "DTHETA": 0.0,
                    "DPSI": 0.0,
                    "S": s_val
                }
                columns = ["NAME", "DX", "DY", "DS", "DPHI", "DTHETA", "DPSI", "S"]
                row_full = {col: 0.0 for col in columns}
                row_full.update({
                    "NAME": name,
                    misalignment_type: float(value_mm),
                    "S": s_val
                })
                self.error_df = pd.concat([self.error_df, pd.DataFrame([row_full])], ignore_index=True)


    def save_to_tfs(self, filename, origin=None):
        from datetime import datetime
        now = datetime.now()
        date_str = now.strftime("%d/%m/%y")
        time_str = now.strftime("%H.%M.%S")
        if origin is None:
            origin = "cpymad"
    
        header_lines = [
            '@ NAME             %06s "EFIELD"',
            '@ TYPE             %06s "EFIELD"',
            '@ TITLE            %08s "no-title"',
            f'@ ORIGIN           %16s "{origin}"',
            f'@ DATE             %08s "{date_str}"',
            f'@ TIME             %08s "{time_str}"',
        ]
    
        # Exact MAD-X column order
        col_names = []
        for i in range(21):
            col_names.append(f'K{i}L')
            col_names.append(f'K{i}SL')
        col_names += [
            'DX', 'DY', 'DS', 'DPHI', 'DTHETA', 'DPSI',
            'MREX', 'MREY', 'MREDX', 'MREDY', 'AREX', 'AREY',
            'MSCALX', 'MSCALY', 'RFM_FREQ', 'RFM_HARMON', 'RFM_LAG'
        ]
        for i in range(21):
            col_names.append(f'P{i}L')
            col_names.append(f'P{i}SL')
    
        # Exact formatting for header and type lines
        col_headers = "* NAME                        " + " ".join(f"{col:<12}" for col in col_names)
        col_types = "$ %s                          " + " ".join("%le".rjust(12) for _ in col_names)
    
        df = self.error_df.copy().sort_values("S", na_position="last")
        if "S" in df.columns:
            df = df.drop(columns=["S"])
    
        with open(filename, "w") as f:
            for line in header_lines:
                f.write(line + "\n")
            f.write(col_headers + "\n")
            f.write(col_types + "\n")
    
            for _, row in df.iterrows():
                name = f'"{row["NAME"].upper()}"'
                line = f' {name:<28}'
    
                for col in col_names:
                    val = row.get(col, 0.0)
                    try:
                        num = float(val)
                    except Exception:
                        num = 0.0
    
                    if col in ['DX', 'DY', 'DS']:
                        num *= 1e-3  # mm → m
                    elif col in ['DPHI', 'DTHETA', 'DPSI']:
                        num *= 1e-3  # mrad → rad
    
                    line += f"{num:.12f} "
    
                f.write(line.rstrip() + "\n")
    
    def add_misalignments_from_dataframe(self, df):
        """
        Process a DataFrame with survey data and apply DY and DPSI misalignments.
    
        Expected columns:
        - magnet, S_start, S_end, S_centre, angle, offset_start, offset_end, offset_centre
    
        Rules:
        - 'Dipole #' → calls add_dipole_misalignment("sp#_dip", ...)
        - 'QD #'     → calls add_quadrupole_misalignment("sp#_qd", ...)
        - 'QF #'     → calls add_quadrupole_misalignment("sp#_qf", ...)
        - 'QC #'     → calls add_quadrupole_misalignment("sp#_qds", ...)
        """
        def map_name_and_type(magnet):
            if magnet.startswith("Dipole "):
                return f"sp{magnet.split()[-1]}_dip", "dipole"
            elif magnet.startswith("QD "):
                return f"sp{magnet.split()[-1]}_qd", "quad"
            elif magnet.startswith("QF "):
                return f"sp{magnet.split()[-1]}_qf", "quad"
            elif magnet.startswith("QC "):
                return f"sp{magnet.split()[-1]}_qds", "quad"
            else:
                raise ValueError(f"Unrecognised magnet label: {magnet}")
    
        df = df.copy()
        df[["name", "type"]] = df["magnet"].apply(lambda m: pd.Series(map_name_and_type(m)))
    
        for _, row in df.iterrows():
            dy = row["offset_centre"]
            dpsi = row["angle"]
            if row["type"] == "dipole":
                self.add_dipole_misalignment(row["name"], "DY", dy)
                self.add_dipole_misalignment(row["name"], "DPSI", dpsi)
            else:
                self.add_quadrupole_misalignment(row["name"], "DY", dy)
                self.add_quadrupole_misalignment(row["name"], "DPSI", dpsi)

class BPMFitResultsLoader:
    def __init__(self, filepath, reverse_co=False):
        self.filepath = filepath
        self.reverse_co = reverse_co
        self.df = self._load_file()
        
    def _load_file(self):
        return pd.read_csv(self.filepath, sep=r"\s+")
    
    def get_bpm_list(self):
        return self.df['bpm'].unique().tolist()

    def get_available_parameters(self):
        return self.df.columns

    def get_parameter(self, param_name, bpm=None):
        if param_name not in self.df.columns:
            raise ValueError(f"Parameter '{param_name}' not found in data.")
        if bpm:
            return self.df[self.df['bpm'] == bpm][param_name].values
        return self.df[param_name].values

    def get_plane_parameters(self, plane="H", parameters=None):
        if parameters is None:
            parameters = self.get_available_parameters()
        missing = [p for p in parameters if p not in self.df.columns]
        if missing:
            raise ValueError(f"Missing parameters in file: {missing}")
        df_plane = self.df[self.df['plane'] == plane].copy()
        return df_plane[["bpm"] + parameters]

    def get_co_data(self, plane="H", twiss_df=None):
        required_columns = ["closed_orbit_mm", "closed_orbit_mm_err"]
        missing = [col for col in required_columns if col not in self.df.columns]
        if missing:
            raise ValueError(f"Missing required CO data columns: {missing}")
        
        df_plane = self.df[self.df['plane'] == plane].copy()
    
        if self.reverse_co:
            df_plane["closed_orbit_mm"] *= -1
    
        # Optional Twiss 's' lookup by partial BPM name
        if twiss_df is not None:
            twiss_lookup = twiss_df.copy()
            twiss_lookup["name_stripped"] = twiss_lookup["name"].astype(str).str.lower()
            bpm_to_s = {}
    
            for bpm in df_plane["bpm"]:
                bpm_key = bpm.lower()
                matched = twiss_lookup[twiss_lookup["name_stripped"].str.contains(bpm_key)]
                if len(matched) == 1:
                    bpm_to_s[bpm] = matched["s"].values[0]
                else:
                    bpm_to_s[bpm] = float("nan")  # ambiguous or not found
    
            df_plane["s"] = df_plane["bpm"].map(bpm_to_s)
    
        return df_plane[["bpm", "closed_orbit_mm", "closed_orbit_mm_err"] + (["s"] if "s" in df_plane.columns else [])]

    def get_dataframe(self):
        return self.df

# %%
def cpymad_set_correctors(madx_instance, cpymad_logfile, corrector_dict, max_E=800., time=0.0):
    """
    Applies the vertical corrector kick values (converted from Amperes to mrad) 
    to a cpymad MAD-X instance.

    Parameters:
    madx_instance (Madx): An instance of cpymad's MAD-X.
    cpymad_logfile (str): Path to the cpymad log file (not used in function).
    corrector_dict (dict): Dictionary with keys as MAD-X variable names 
                           and values as the programmed kicks in Amperes.
    """

    for key, amps in corrector_dict.items():
        # Extract plane ('V' from 'vd1') and super-period (from 'rX' where X is 0, 2, etc.)
        sp = int(key[1])  # Extract the second character as the super-period
        plane = 'V' if 'vd' in key else 'H'  # Determine plane from key name

        # Convert kick from Amperes to milliradians
        kick_mrad = calculate_corrector_kick(amps, max_E, time, plane, sp)

        # Print key, amps, and converted kick
        print(f"{key}: {amps:.6f} A -> {kick_mrad:.6f} mrad")

        # Apply the converted kick to MAD-X
        kick_mrad *= 1E-3 # convert from millirad to radians
        madx_instance.input(f"{key} := {kick_mrad};")


# %%
def cpymad_apply_and_check_error_table(madx_instance, error_file, original_df, atol=1e-10, rtol=1e-12):
    """
    Apply the error table to MAD-X, extract the resulting table, and compare with the original error_df.

    Parameters:
        madx_instance : cpymad.madx.Madx
            Active MAD-X instance.
        error_file : str
            Path to the .tfs error table file.
        original_df : pd.DataFrame
            The builder.error_df to compare against (values in mm/mrad).
        atol : float
            Absolute tolerance for comparison (default: 1e-10).
        rtol : float
            Relative tolerance for comparison (default: 1e-12).

    Returns:
        bool
            True if the tables match within tolerance, otherwise raises AssertionError.
    """
    # Apply error table
    cpymad_apply_error_table(madx_instance, error_file)

    # Get MAD-X table in m/rad
    madx_df = get_madx_table_df(madx_instance, nonzero=False)
    madx_df = madx_df[["name", "dx", "dy", "ds", "dphi", "dtheta", "dpsi"]].copy()
    madx_df["name"] = madx_df["name"].str.lower()

    # Scale original_df from mm/mrad → m/rad
    original = original_df[["NAME", "DX", "DY", "DS", "DPHI", "DTHETA", "DPSI"]].copy()
    original.columns = [c.lower() for c in original.columns]
    original["name"] = original["name"].str.lower()

    for col in ["dx", "dy", "ds", "dphi", "dtheta", "dpsi"]:
        original[col] = original[col].astype(float) * 1e-3

    # Sort both DataFrames
    madx_sorted = madx_df.sort_values("name").reset_index(drop=True)
    original_sorted = original.sort_values("name").reset_index(drop=True)

    # Compare with tolerance
    pd.testing.assert_frame_equal(madx_sorted, original_sorted, check_dtype=False, rtol=rtol, atol=atol)
    return True

def cpymad_apply_error_table(madx_instance, error_table_file):
    """
    Apply a MAD-X error table file using madx_instance.input().

    Parameters:
        madx_instance: The cpymad.Madx instance.
        error_table_file: Path to the .tfs error table file.
    """
    madx_instance.input(f'READMYTABLE, file="{error_table_file}", table=efield;')
    madx_instance.input('SETERR, TABLE=efield;')

def get_madx_table_df(madx, table_name="efield", nonzero=True):
    """
    Extracts a MAD-X table as a pandas DataFrame.

    Parameters:
        madx : cpymad.madx.Madx
            The active MAD-X instance.
        table_name : str
            Name of the MAD-X table to extract (default: "efield").
        nonzero : bool
            If True, return only rows where any of the key misalignment/rotation columns are non-zero.

    Returns:
        pd.DataFrame
            DataFrame containing the filtered or full table data.
    """
    if table_name not in list(madx.table):
        raise ValueError(f"MAD-X table '{table_name}' not found. Available tables: {list(madx.table)}")

    raw_table = getattr(madx.table, table_name)
    raw_data = raw_table.copy()
    df = pd.DataFrame(raw_data, columns=raw_data.keys())

    if nonzero:
        cols = ["ds", "dx", "dy", "dtheta", "dphi", "dpsi"]
        df = df[(df[cols] != 0).any(axis=1)][["name"] + cols]

    return df
def df_to_correction_dict(df, plane='v'):
    correction_dict = {}
    plane = plane.lower()

    for _, row in df.iterrows():
        name = row['NAME'].lower()
        if plane == 'v' and 'vd' in name:
            value = row['PY.CORRECTION'] * 1E3
        elif plane == 'h' and 'hd' in name:
            value = row['PX.CORRECTION'] * 1E3
        else:
            continue

        # Extract something like 'r0vd1' from 'sp0_r0vd1'
        parts = name.split('_')
        if len(parts) > 1:
            shortname = parts[1]
        else:
            shortname = name

        correction_dict[f"{shortname}_kick"] = round(value, 3)

    return correction_dict

def calculate_corrector_kick(current_Amps, max_E, time, plane='H', sp=0):
    """
    Returns the corrector steering kick in milliradians given the desired current in amperes.

    Parameters:
    current_Amps (float): Desired current in Amperes.
    max_E (float): Maximum energy.
    time (float): Measurement time.
    plane (str): 'H' for horizontal or 'V' for vertical.
    sp (int): Super-period number.

    Returns:
    float: Corrector current in amperes.
    """
    sp_list = [0, 2, 3, 4, 5, 7, 9]
    if sp not in sp_list:
        print('calculate_corrector_current:: selected super-period has no steering magnet')
        exit(0)

    # Calibration provided by HVC 30.09.22
    calibration_data = {
        '0H': 0.08350, '2H': 0.09121, '3H': 0.08, '4H': 0.06600,
        '5H': 0.07780, '7H': 0.07580, '9H': 0.07660, '0V': 0.04620,
        '2V': 0.04330, '3V': 0.05210, '4V': 0.04770, '5V': 0.05400,
        '7V': 0.05220, '9V': 0.04510
    }

    df = synchrotron_energy_data(max_E, time)

    h_list = ['h', 'H', 'horizontal', 'Horizontal']
    key = f"{sp}{'H' if plane in h_list else 'V'}"
    
    # Compute the kick in milliradians
    kick_mrad  = current_Amps /( df['Rigidity [Tm]'].iloc[0] / calibration_data[key])

    return round_sig(kick_mrad,7)


# %%
def calculate_corrector_current(kick_mrad, max_E, time, plane='H', sp=0):
    """
    Returns the corrector current in amperes given the desired steering kick in milliradians.

    Parameters:
    kick_mrad (float): Desired steering kick in milliradians.
    max_E (float): Maximum energy.
    time (float): Measurement time.
    plane (str): 'H' for horizontal or 'V' for vertical.
    sp (int): Super-period number.

    Returns:
    float: Corrector current in amperes.
    """
    sp_list = [0, 2, 3, 4, 5, 7, 9]
    if sp not in sp_list:
        print('calculate_corrector_current:: selected super-period has no steering magnet')
        exit(0)

    # Calibration provided by HVC 30.09.22
    calibration_data = {
        '0H': 0.08350, '2H': 0.09121, '3H': 0.08, '4H': 0.06600,
        '5H': 0.07780, '7H': 0.07580, '9H': 0.07660, '0V': 0.04620,
        '2V': 0.04330, '3V': 0.05210, '4V': 0.04770, '5V': 0.05400,
        '7V': 0.05220, '9V': 0.04510
    }

    df = synchrotron_energy_data(max_E, time)

    h_list = ['h', 'H', 'horizontal', 'Horizontal']
    key = f"{sp}{'H' if plane in h_list else 'V'}"

    # Compute the current in amperes
    amps = kick_mrad * df['Rigidity [Tm]'].iloc[0] / calibration_data[key]

    return round_sig(amps,7)


# %%
def convert_kicks_to_currents(kicks_dict, max_E=800., time=0.0):
    """
    Converts a dictionary of corrector kicks in milliradians to corrector settings in amperes.

    Parameters:
    kicks_dict (dict): Dictionary where keys are corrector names and values are kicks in milliradians.
    max_E (float): Maximum energy, default is 800.
    time (float): Measurement time, default is 0.0.

    Returns:
    dict: A new dictionary with the same keys but values converted to amperes.
    """
    current_dict = {}

    for key, kick_mrad in kicks_dict.items():
        # Extract plane ('V' from 'vd1') and super-period (from 'rX' where X is 0, 2, etc.)
        sp = int(key[1])  # Extract super-period from the second character
        plane = 'V' if 'vd' in key else 'H'  # Determine plane based on key name

        # Convert the kick from milliradians to amperes
        amps = calculate_corrector_current(kick_mrad, max_E, time, plane, sp)

        # Store the converted value in the new dictionary
        current_dict[key] = amps

    return current_dict


# %%
def convert_currents_to_kicks(currents_dict, max_E=800., time=0.0):
    """
    Converts a dictionary of corrector settings in amperes to corrector kicks in milliradians.

    Parameters:
    currents_dict (dict): Dictionary where keys are corrector names and values are kicks in milliradians.
    max_E (float): Maximum energy, default is 800.
    time (float): Measurement time, default is 0.0.

    Returns:
    dict: A new dictionary with the same keys but values converted to amperes.
    """
    current_dict = {}

    for key, current_Amp in currents_dict.items():
        # Extract plane ('V' from 'vd1') and super-period (from 'rX' where X is 0, 2, etc.)
        sp = int(key[1])  # Extract super-period from the second character
        plane = 'V' if 'vd' in key else 'H'  # Determine plane based on key name

        # Convert the kick from milliradians to amperes
        amps = calculate_corrector_kick(current_Amp, max_E, time, plane, sp)

        # Store the converted value in the new dictionary
        current_dict[key] = amps

    return current_dict


def cpymad_set_isis_cycle_time(madx_instance, max_E, time):
    # Ensure time is a float and in valid increments
    if not isinstance(time, float) or time < 0.0 or time > 10.0 or (time * 10) % 5 != 0:
        print(f"Error: time must be a float between 0.0 and 10.0 in 0.5 increments. Received: {time}")
        return time

    # Generate dataframe of synchrotron energy and related info
    energy_df = synchrotron_energy_df(max_E, intervals=20)

    # store some values for this time point
    try:
        energy = energy_df[energy_df['Time [ms]'] == time]['Energy [MeV]'].iloc[0]
        pc = energy_df[energy_df['Time [ms]'] == time]['Momentum [GeV/c]'].iloc[0]
    except IndexError:
        print(f"Error: No matching time value found in energy dataframe for time = {time} ms")
        return time

    # set the beam to this energy in cpymad
    madx_instance.input(f'beam, particle = proton, pc = {pc};')

    # print confirmation
    print(f'ISIS cpymad run, energy set to {energy} MeV, pc = {pc}')
    return time

# Streamlit title and description

def apply_correctors(madx, twiss_input, corrector_dict, max_E, cycle_time):
    cpymad_set_correctors(madx, cpymad_logfile, corrector_dict, max_E, cycle_time)
    return cpymad_madx_twiss(madx, cpymad_logfile, sequence_name)


def apply_misalignments(madx, twiss_input, misalignments_file):
    df_misalignments = pd.read_csv(misalignments_file, sep="\t")

    error_table_builder = cpymad_ErrorTableBuilder(twiss_input)
    error_table_builder.add_misalignments_from_dataframe(df_misalignments)
    error_file = "uploaded_vertical_misalignments.tfs"
    error_table_builder.save_to_tfs(error_file, origin="UPLOADED")

    cpymad_apply_and_check_error_table(madx, error_file, error_table_builder.error_df)
    return cpymad_madx_twiss_nocheck(madx, cpymad_logfile, sequence_name)

def generate_orbit_plot(twiss_data, epsilon, title_suffix="", overlay_data=None, xlimits=None):
    qx = madx.table.summ.q1[0]
    qy = madx.table.summ.q2[0]
    plot_title = f"{sequence_name} Q1={qx:.3f}, Q2={qy:.3f} {title_suffix}"

    fig = make_subplots(
        rows=2, cols=1,
        shared_xaxes=True,
        vertical_spacing=0.02,
        row_heights=[0.25, 0.25],
        specs=[[{"type": "scatter"}],
               [{"type": "scatter"}]]
    )

    # Calculating envelopes

    betx_array = twiss_data["betx"]
    bety_array = twiss_data["bety"]

    sigma_x = np.sqrt(betx_array*epsilon)
    sigma_y = np.sqrt(bety_array*epsilon)

    env_plus_x = twiss_data['x']  + sigma_x
    env_minus_x = twiss_data['x'] - sigma_x

    env_plus_y = twiss_data['y'] + sigma_y
    env_minus_y = twiss_data['y'] - sigma_y

    # Row 1: horizontal orbit
    fig.add_trace(go.Scatter(x=twiss_data['s'], y=twiss_data['x'] *1e3, mode='lines', name='Horizontal Closed Orbit', line=dict(color='black')), row=1, col=1)
    # Row 1: Envelope
    fig.add_trace(go.Scatter(x=twiss_data['s'], y=env_plus_x*1e3, mode='lines', name='Envelope', line=dict(color='orange')), row=1, col=1)
    fig.add_trace(go.Scatter(x=twiss_data['s'], y=env_minus_x*1e3, mode='lines', line=dict(color='orange')), row=1, col=1)
    # Row 2: vertical orbit (with optional overlay)
    fig.add_trace(go.Scatter(x=twiss_data['s'], y=twiss_data['y'] *1e3, mode='lines', name='Vertical Closed Orbit', line=dict(color='black')), row=2, col=1)
    # Row 2: Envelope
    fig.add_trace(go.Scatter(x=twiss_data['s'], y=env_plus_y*1e3, mode='lines', line=dict(color='orange')), row=2, col=1)
    fig.add_trace(go.Scatter(x=twiss_data['s'], y=env_minus_y*1e3, mode='lines', line=dict(color='orange')), row=2, col=1)


    if overlay_data is not None:
        fig.add_trace(go.Scatter(
            x=overlay_data['s'],
            y=overlay_data['closed_orbit_mm'],
            mode='markers',
            name='Measured vertical orbit',
            marker=dict(size=8, color='red')
        ), row=2, col=1)

    # Layout
    fig.update_layout(title=plot_title, height=800, width=1000,
                      legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
                      template="plotly_white")

    fig.update_yaxes(title_text='x [mm]', row=1, col=1)
    fig.update_yaxes(title_text='y [mm]', row=2, col=1)
    fig.update_xaxes(title_text='s [m]', row=2, col=1)
    if xlimits is not None:
        fig.update_xaxes(range=xlimits)

    st.plotly_chart(fig)

#Loading twiss dataset, initialisation of cycle time and energy
plot_folder = 'Orbit_Correction_Plots'
make_directory(plot_folder)
st.title("BPM Simulation for Accelerator Physics")
st.write("""
    This application allows you to simulate and visualize the beam position monitor (BPM) results in an accelerator. 
    You can upload your BPM data, adjust corrector settings, and visualize the orbit data.
""")
max_E = 800 # 800 MeV
cycle_time = 0.0 
time_point = cpymad_set_isis_cycle_time(madx, max_E, cycle_time)
twiss_0 = cpymad_madx_twiss(madx, cpymad_logfile, sequence_name)

st.write("Twiss Data Preview:", twiss_0.head())


##Logo 
import base64
with open("ukri-stfc-square-logo.png", "rb") as f:
    data = base64.b64encode(f.read()).decode("utf-8")

    st.sidebar.markdown(
        f"""
        <div style="display:table;margin-top:-20%;margin-left:20%;">
            <img src="data:image/png;base64,{data}" width="150" height="150">
        </div>
        """,
        unsafe_allow_html=True)

# Input fields for corrector settings
st.sidebar.header("Corrector Settings")

# Corrector input: Default current values
corrector_values = {
    "r0vd1_kick": 16.3,
    "r2vd1_kick": -33.0,
    "r3vd1_kick": 30,
    "r4vd1_kick": 22.8,
    "r5vd1_kick": -21.9,
    "r7vd1_kick": -38.1,
    "r9vd1_kick": 37.8
}

corrector_currents = {}

# FOR VD: Convert the kicks to currents if necessary
corrector_currents_converted = convert_kicks_to_currents(corrector_currents, max_E, cycle_time)
# Run MAD-X and get the closed orbit with the new corrector settings
v_corrected_dict = {
    "r0vd1_kick": -0.073474149,
    "r2vd1_kick": 0.307515168,
    "r3vd1_kick": -0.257864760,
    "r4vd1_kick": -0.019923798,
    "r5vd1_kick": -0.319146002,
    "r7vd1_kick": 0.549252123,
    "r9vd1_kick": 0.603084818,
}
v_corrector_currents = convert_kicks_to_currents(v_corrected_dict)
v_corrector_kicks = convert_currents_to_kicks(v_corrector_currents)
v_corrector_currents_minus_0p4ms = {
 'r0vd1_kick': 16.3,
 'r2vd1_kick': -33.0,
 'r3vd1_kick': 30,
 'r4vd1_kick': 22.8,
 'r5vd1_kick': -21.9,
 'r7vd1_kick': -38.1,
 'r9vd1_kick': 37.8}


# FOR HD
h_corrector_currents_minus_0p4ms = {
 'r0hd1_kick': 5,
 'r2hd1_kick': -5,
 'r3hd1_kick': 10,
 'r4hd1_kick': 2,
 'r5hd1_kick': -2,
 'r7hd1_kick': -22,
 'r9hd1_kick': 15}


#TODO: Remove all corrections
#TODO: Error table: choose file to edit 

# # Sidebar checkboxes
#apply_corr = st.sidebar.checkbox("Apply Vertical Correctors")
# apply_mis = st.sidebar.checkbox("Apply Misalignments")

#User can enter their own beam emittance
epsilonOrig= float(st.sidebar.number_input("Enter Beam emittance in mmmrad (default is 300 )", 
                                       min_value = float(50), 
                                       max_value = float(500),
                                       value = float(300)))
if epsilonOrig ==0:
    epsilon = 300 * 1e-6
else:
    epsilon = epsilonOrig * 1e-6
st.sidebar.write("The beam emittance is ", str(epsilonOrig), "mmmrad")

# dipole correctors (change orbit)
apply_hd = st.sidebar.checkbox("Apply Horizontal Dipole") ##x
apply_vd = st.sidebar.checkbox("Apply Vertical Dipole")  ##y  ### <--- TO BE CONTINUED
# trim quads (change tune - focusing)
apply_tunes = st.sidebar.checkbox("Apply Tunes")
apply_harmonics = st.sidebar.checkbox("Apply Harmonics") 

uploaded_mis_file = None
# if apply_mis:
#     uploaded_mis_file = st.sidebar.file_uploader("Upload misalignments file (.txt)", type=["txt"])

uploaded_bpm_file = st.sidebar.file_uploader("Upload BPM fit results (.txt)", type=["txt"])

# Simulation Controls
max_E = 800
cycle_time_slider = st.sidebar.slider("Cycle Time (ms)", min_value=0.0, max_value=10.0, step=0.5, value=0.0)

# Optional xlimit restriction
restrict_xaxis = st.sidebar.checkbox("Restrict x-axis domain to 1 super-period")
xlimits = [16.336282 * 4, 16.336282 * 5] if restrict_xaxis else None

twiss_current = twiss_0.copy()

## Vertical corrections (not just a straight line)
if apply_vd:
    twiss_current = apply_correctors(madx, twiss_current, v_corrector_currents_minus_0p4ms, max_E, cycle_time_slider)

if apply_hd:
    twiss_current = apply_correctors(madx, twiss_current, h_corrector_currents_minus_0p4ms, max_E, cycle_time_slider)

if apply_tunes:
    df = getValues()
    Qh = df['time'==time_point]['x']
    Qv = df['time'==time_point]['y']
    set_tune_DW(madx, cpymad_logfile, Qh, Qv, time_point)

if apply_harmonics:
    madx.globals['D7COS'] = get_harmonic("D7COS", time_point)
    madx.globals['D8COS'] = get_harmonic("D8COS", time_point)
    madx.globals['F8COS'] = get_harmonic("F8COS", time_point)
    madx.globals['D7SIN'] = get_harmonic("D7SIN", time_point)
    madx.globals['D8SIN'] = get_harmonic("D8SIN", time_point)
    madx.globals['F8SIN'] = get_harmonic("F8SIN", time_point)



# if apply_mis and uploaded_mis_file:
#     twiss_current = apply_misalignments(madx, twiss_current, uploaded_mis_file)
#     v_corrector_currents_off = {
#     'r0vd1_kick': 0,
#     'r2vd1_kick': 0,
#     'r3vd1_kick': 0,
#     'r4vd1_kick': 0,
#     'r5vd1_kick': 0,
#     'r7vd1_kick': 0,
#     'r9vd1_kick': 0}

#     #calculate IBO corrections 
#     cpymad_set_correctors(madx, cpymad_logfile, v_corrector_currents_off, max_E, cycle_time)
#     madx.command.usekick(sequence=sequence_name, status="on", pattern="^R.*VD.*")
#     madx.command.usemonitor(sequence=sequence_name, status="on", class_="monitor")
#     bare_madx_twiss_file = sequence_name +'_madx_twiss_bare.tfs'
#     madx.input('set, format="12.12f"')
#     madx.input('select, flag=twiss, column=keyword, name, s, l, betx, alfx, mux, bety, alfy, muy, x, px, y, py, t, pt, dx, dpx, dy, dpy, wx, phix, dmux, wy, phiy, dmuy, ddx, ddpx, ddy, ddpy, r11, r12, r21, r22, energy, angle, k0l, k0sl, k1l, k1sl, k2l, k2sl, k3l, k3sl, k4l, k4sl, k5l, k5sl, k6l, k6sl, k7l, k7sl, k8l, k8sl, k9l, k9sl, k10l, k10sl, ksi, hkick, vkick, tilt, e1, e2, h1, h2, hgap, fint, fintx, volt, lag, freq, harmon, slot_id, assembly_id, mech_sep, kmax, kmin, calib, polarity, alfa, beta11, beta12, beta13, beta21, beta22, beta23, beta31, beta32, beta33, alfa11, alfa12, alfa13, alfa21, alfa22, disp1, disp2, disp3, disp4')
#     twiss_bare = madx.twiss(sequence=sequence_name, file=bare_madx_twiss_file, table='bare').dframe()
#     madx.command.correct(model='bare',sequence=sequence_name, plane="y", flag="ring", error=1e-7, mode='svd', cond=1, corzero=1, monerror=0, monscale=0, clist=clist_file_y, mlist=mlist_file_y)   
#     corrected_madx_twiss_file = sequence_name +'_madx_twiss_corrected_.tfs'
#     madx.input('set, format="12.12f"')
#     madx.input('select, flag=twiss, column=keyword, name, s, l, betx, alfx, mux, bety, alfy, muy, x, px, y, py, t, pt, dx, dpx, dy, dpy, wx, phix, dmux, wy, phiy, dmuy, ddx, ddpx, ddy, ddpy, r11, r12, r21, r22, energy, angle, k0l, k0sl, k1l, k1sl, k2l, k2sl, k3l, k3sl, k4l, k4sl, k5l, k5sl, k6l, k6sl, k7l, k7sl, k8l, k8sl, k9l, k9sl, k10l, k10sl, ksi, hkick, vkick, tilt, e1, e2, h1, h2, hgap, fint, fintx, volt, lag, freq, harmon, slot_id, assembly_id, mech_sep, kmax, kmin, calib, polarity, alfa, beta11, beta12, beta13, beta21, beta22, beta23, beta31, beta32, beta33, alfa11, alfa12, alfa13, alfa21, alfa22, disp1, disp2, disp3, disp4')
#     twiss_corrected = madx.twiss(sequence=sequence_name, file=corrected_madx_twiss_file, table='corrected').dframe()
    
#     #plotting of IBO vs corrections
#     fig4 = go.Figure()

#     fig4.add_trace(go.Scatter(
#         x=twiss_current.s,
#         y=twiss_current.y * 1e3,
#         mode='lines',
#         name='Inferred bare orbit (IBO) from misalignments'
#     ))

#     fig4.add_trace(go.Scatter(
#         x=twiss_corrected.s,
#         y=twiss_corrected.y * 1e3,
#         mode='lines',
#         name='IBO + MAD-X Correction'
#     ))

#     fig4.update_layout(
#         title='Comparison of bare orbit inferred from survey misalignments, and corrected orbit',
#         xaxis_title='S [m]',
#         yaxis_title='y [mm]',
#         legend=dict(title=None),
#         template='simple_white',
#         margin=dict(t=60, b=40, l=60, r=20),
#     )
#     fig4.update_xaxes(showgrid=True, gridwidth=0.5, gridcolor='lightgrey')
#     fig4.update_yaxes(showgrid=True, gridwidth=0.5, gridcolor='lightgrey')

#     st.plotly_chart(fig4)

# # Load overlay BPM if provided
# bpm_overlay = None
# if uploaded_bpm_file:
#     bpm_bare = BPMFitResultsLoader(uploaded_bpm_file, reverse_co=False)
#     bpm_overlay = bpm_bare.get_co_data('V', twiss_0)

# Generate the main plot
# generate_orbit_plot(twiss_current, title_suffix="", overlay_data=bpm_overlay, xlimits=xlimits)
generate_orbit_plot(twiss_current, epsilon, title_suffix="", xlimits=xlimits)