import pandas as pd

def get_EPICS_Tune(cycletime:str):
    time_periods = [-.4, -.2, 0, .5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8.1, 9, 9.8, 10.1]
    qx_array = [4.315, 4.270, 4.270, 4.250, 4.235, 4.205, 4.170, 4.190, 4.18, 4.18, 4.18, 4.17, 4.165, 4.165, 4.165, 4.18, 4.18, 4.175, 4.170, 4.190, 4.18, 4.18]
    qy_array = [3.82, 3.82, 3.81, 3.805, 3.800, 3.825, 3.680, 3.680, 3.69, 3.7, 3.7, 3.695, 3.695, 3.695, 3.692, 3.69, 3.680, 3.680, 3.69, 3.7, 3.82, 3.82]

    data = {
        "PV":["", ""],
        "Q_request":[],
        "Last_change":["", ""]
    }

    for i in range(0, len(time_periods)):
        if cycletime == time_periods[i]:
            data["Q_request"].append(qx_array[i])
            data["Q_request"].append(qy_array[i])
    
    df = pd.DataFrame.from_dict(data)
    return df

def get_EPICS_Beta(twiss, cycletime:str):
    time_periods = [-.4, -.2, 0, .5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8.1, 9, 9.8, 10.1]
    betx_array = twiss["betx"]
    bety_array = twiss["bety"]

    data = {
        "PV":["", ""],
        "Q_request":[],
        "Last_change":["", ""]
    }

    for i in range(0, len(time_periods)):
        if cycletime == time_periods[i]:
            data["Q_request"].append(betx_array[i])
            data["Q_request"].append(bety_array[i])
    
    df = pd.DataFrame.from_dict(data)
    return df
    


# Takes a cycle time
# PV - DWTRIM::H_Q:AT_TIME:XMS
#      DWTRIM::V_Q:AT_TIME:XMS
