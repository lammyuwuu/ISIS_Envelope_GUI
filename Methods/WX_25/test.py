import pandas as pd
import os



def get_harmonic(pv_name: str, cycletime: str | float) -> float:
    """Return rows matching a specific CycleTime and PV value."""
    # Get the directory of the current script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    csv_path = os.path.join(script_dir, "..", "Collected_EPICS_data", "get_EPICS_Harmonics_full_cycle.dat")

    harmonic_data = pd.read_csv(csv_path)

    harmonic_index = harmonic_data["PV"].str.contains(f"DWTRIM::{pv_name}:AT_TIME:{cycletime}MS")
    return harmonic_data[harmonic_index]["Harmonic"].to_list()[0]

# Example usage:
result = get_harmonic("D7COS", "0")
print(result)
