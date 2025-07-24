import time
import numpy as np
import pandas as pd
from pandas import DataFrame # Used for TypeHinting
import pandas as pd
from typing import List # Used for TypeHinting
import datetime
import requests
import re 
import os
from math import log10, floor

'''A collection of functions originally made by the 2024 Work Experience group. Has been edited and formatted to work 1 year later.
Note: get_EPICS_**** functions are dependent on p4p to access live data, which is dependent on C++ compilers, which are only installed on 2 of the laptops.
Set the environment variable os.environ["ACCESS_LIVE_DATA"] = 'False' before importing to avoid importing p4p.'''

def get_historical_data(pv_name: str, start_time: datetime.datetime, end_time: datetime.datetime, archiver_addr: str = "http://athena.isis.rl.ac.uk:9505") -> DataFrame:
    """Function to get Program Variables from the archive between given dates."""
    #Convert DataTime object to 
    start_time = start_time.strftime("%Y-%m-%dT%H:%M:%S.%fZ")
    end_time = end_time.strftime("%Y-%m-%dT%H:%M:%S.%fZ")

    
    #API String to get data from a set PV and its time window of intrest
    data_str_req = f'{archiver_addr}/data?pv={pv_name}&from={start_time}&to={end_time}'
    print(data_str_req)
    data_meta = requests.get(data_str_req).json()

    # Convert the list of dictionaries to a DataFrame
    df = pd.DataFrame(data_meta)

    # Combine 'secs' and 'nanos' to create a datetime index
    df.insert(0,'timestamp', pd.to_datetime(df['secs'], unit='s') + pd.to_timedelta(df['nanos'], unit='ns'))

    # Drop the 'secs' and 'nanos' columns if they are no longer needed
    df.drop(columns=['secs', 'nanos'], inplace=True)
    
    return df

def date_to_unix(date: datetime.datetime) -> float:
    """Helper function to reformat dates."""
    return date.timestamp()

def get_archive_trim_quads(t_cycle: str = "5", magnet_type: str = "QTD", year: str = "2024", month: str = "10", day: str = "05", t_start: str = "01:01:02", t_end: str = "23:59:59", archiver_addr: str = "http://athena.isis.rl.ac.uk:9505") -> DataFrame:
    """Function to get the current at the trim quads from the archive, given the cycle, magnet type, and date. Calls get_historical_data()."""
    current_time = datetime.datetime.now() #loads the current time and an arbitrarily set start time which is crucial for later on
    start_time = datetime.datetime(2022, 6, 20) # when we will be loading up data
    time_periods = ["-.6","-.5","-.4","-.3","-.2","-.1","0",".5","1","1.5","2","2.5","3","3.5","4","4.5","5","5.5","6","7","8","9","10"]
    regex = "^(?:[01][0-9]|2[0-3]):[0-5][0-9](?::[0-5][0-9])?$" # this is simply some data validation code

    if str(t_cycle) not in time_periods:
        print("Please enter the t_cycle variable in the format \"-.6\" as a string")
        print("Please enter a valid t_cycle value, please see the values listed below") #make sure you input in a valid time
        print(" | ".join(time_periods))
    t_cycle = str(t_cycle)# allows you to input in the time as an int on certain occasions

    if re.search(t_start,regex): #just some code which ensures the time inputs are in the correct format
        print("The t_start variable is not valid please enter in a value in the format of HH:MM:SS as a string")
    else:
        pass
    
    requested_array = []
    start_timestamp = year+"-"+month+"-"+day+" "+t_start+".000000"#just converting the inputs to datetime format
    end_timestamp = year+"-"+month+"-"+day+" "+t_end+".999999"

    start_timestamp = date_to_unix(datetime.datetime.strptime(start_timestamp,'%Y-%m-%d %H:%M:%S.%f'))
    end_timestamp =  date_to_unix(datetime.datetime.strptime(end_timestamp,'%Y-%m-%d %H:%M:%S.%f')) # converting the datetime format
    # to unix format
    channel_list = []

    count = 0 # loads the valid time inputs below the t_end into a large array 
    full_dataframe = []
    for j in range(10):
        current_channel_string = "DWQ_TEST::R"+str(j)+"QTF:CURRENT:"+t_cycle+"MS"
        channel_list.append(get_historical_data(current_channel_string,start_time,current_time, archiver_addr=archiver_addr))

    count = 0 # loads the valid time inputs below the t_end into a large array 
    full_dataframe = []
    
    for i in range(10):
        for index,row in channel_list[i].iterrows():
            if row.iloc[1] != 0 and date_to_unix(row.iloc[0]) <= end_timestamp:
                full_dataframe.append([date_to_unix(row.iloc[0]),float(row.iloc[1]),"R"+str(i)+"QTD"])
            if date_to_unix(row.iloc[0]) > end_timestamp:
                break


    full_dataframe.sort()#ensures the currents are sorted by time as data is only recorded when the quadripoles are changed
    full_dataframe = full_dataframe[::-1]

    conditions_satisfied_D = [-1]*10
    conditions_satisfied_F = [-1]*10

    # gets the valid times so that each trim quadripole has one value and that value is the most recent one
    valid = True
    for item in full_dataframe:
        if valid == False:
            break

        if item[2][-1] == "D":
            if conditions_satisfied_D[int(item[2][1]) - 1] == -1:
                requested_array.append(item)
                conditions_satisfied_D[int(item[2][1]) - 1] = 1
            if sum(conditions_satisfied_D)+sum(conditions_satisfied_F) == 20:
                valid = False
        else:
            if conditions_satisfied_F[int(item[2][1]) - 1] == -1:
                requested_array.append(item)
                conditions_satisfied_F[int(item[2][1]) - 1] = 1
            if sum(conditions_satisfied_D)+sum(conditions_satisfied_F) == 20:
                valid = False

    for i in range(len(requested_array)):
        requested_array[i] = [requested_array[i][1], (datetime.datetime.fromtimestamp(requested_array[i][0])).replace(second=0, microsecond=0) , t_cycle, requested_array[i][-1]]
    result_dataframe = pd.DataFrame(requested_array,columns=["Current","Last Change","Cycle Time","Trim Quad"])
    result_dataframe = result_dataframe[['Cycle Time', 'Trim Quad', 'Current', 'Last Change']]

    return result_dataframe 

def date_to_unix(date: datetime.datetime | str) -> float:
    """Helper function to reformat date."""
    if isinstance(date, str):
        date = datetime.datetime.strptime(date, '%Y-%m-%d %H:%M:%S.%f')
    return date.timestamp()

def get_archive_correctors(cycleTime: str, axis: str, timeStart: datetime.datetime, timeEnd: datetime.datetime, archiver_addr: str = "http://athena.isis.rl.ac.uk:9505") -> DataFrame: # returns a pandas datatype with headers timestamp, val for current, magnet name and cycle time
    """Function to get corrector values from the archive, given time in cycle, and date. Calls get_historical_data()."""
    returnData = pd.DataFrame()
    for i in range(0,10):
        if i == 1 or i == 8 or i == 6:
            continue
        else:
            print(f"DW{axis}ST_TEST::R{i}{axis}D1:CURRENT:{cycleTime}MS")
            HistoryData = (get_historical_data(f"DW{axis}ST_TEST::R{i}{axis}D1:CURRENT:{cycleTime}MS",timeStart,timeEnd, archiver_addr=archiver_addr))
            
            while HistoryData.empty:
                HistoryData = (get_historical_data(f"DW{axis}ST_TEST::R{i}{axis}D1:CURRENT:{cycleTime}MS",timeStart-datetime.timedelta(month=1),timeEnd, archiver_addr=archiver_addr))
            HistoryData.insert(1,"Magnet Name",f"R{i}{axis}D1")
            returnData = pd.concat([returnData,HistoryData])
    returnData.insert(2,"Cycle Time",cycleTime)
    return returnData

def convert_to_df(pv_data_list: list[list]) -> DataFrame:
    """Helper function to reformat data."""
    column_count = len(pv_data_list[0])
    column_names = []
    for i in range(column_count):
        column_names.append("DataPoint " + str(i))
    df = pd.DataFrame(pv_data_list, columns=column_names)
    return df

def save_to_csv(pv_df: DataFrame, filename: str) -> None: 
    """Helper function to save DataFrame to a csv file."""
    pv_df.to_csv(filename,index=False)
    
def search_pvs(query: str, pvs_addr: str = "http://athena.isis.rl.ac.uk:9505/getPVStatus") -> list[dict]:
    """
    Search for a specified query, using regex, and return a list of PVs

    For example, we say 'TGT1*' we will get a list of PVs that start with 'TGT1'
    """

    url = pvs_addr
    resp = requests.get(url, params={"pv": query, "reporttype": "short"})
    return resp.json()

def get_pv_names(pv_list: list[dict]) -> list[str]:
    """
    Extract names only from the list of PVs returned by the search_pvs function
    """
    pv_name_list = []
    pv_list_len = (len(pv_list))
    for i in range(pv_list_len):
        pv_name_list.append(pv_list[i]["pvName"])
    return pv_name_list

def round_sig(x: float, sig: int = 3) -> float | int:
    """Helper function to round the signal."""
    if x == 0.0: 
        return 0.0
    else: 
        return round(x, sig - int(floor(log10(abs(x)))) - 1)

if os.environ.get("ACCESS_LIVE_DATA"):
    if os.environ["ACCESS_LIVE_DATA"] == "True":
        from p4p.client.thread import Context 
    elif os.environ.get("SUPPRESS_P4P_WARNING"): 
        if os.environ["SUPPRESS_P4P_WARNING"] == "False":
            print("Warning: LIVE DATA UNAVAILABLE. 'get_EPICS_****' functions are dependent on p4p which can only be run on the specific p4p laptops. Will give an error if accessed. Use os.environ['SUPPRESS_P4P_WARNING'] = 'True' to suppress this warning.")
    else:
        print("Warning: LIVE DATA UNAVAILABLE. 'get_EPICS_****' functions are dependent on p4p which can only be run on the specific p4p laptops. Will give an error if accessed. Use os.environ['SUPPRESS_P4P_WARNING'] = 'True' to suppress this warning.")

else:
    from p4p.client.thread import Context
import datetime

def get_EPICS_Horizontal_Correctors(cycletime: str, filename: str | None = None, save: bool = False, pvs_addr: str = "http://athena.isis.rl.ac.uk:9505/getPVStatus", live_addr: str = "infra.isis.rl.ac.uk:7075", max_tries: int = 3) -> DataFrame: #in MS
    """Function to get live Horizontal Corrector data from EPICS. Returns a DataFrame, given the time in cycle. Will save to filename if save is True. If filename is None, will automatically set the filename.  This function is dependent on p4p, see module description. This function will not work unless p4p imported."""
    os.environ["EPICS_PVA_NAME_SERVERS"] = live_addr
    ctxt = Context("pva")
    if filename is None: filename='EPICS_Tune_HC_'+str(cycletime)+'.dat'
        
    time_periods = ["-.4", "-.2", "0", ".5", "1", "1.5", "2", "2.5", "3", "3.5", "4", "4.5", "5", "5.5", "6", "6.5", "7", "7.5", "8", "8.5", "9", "10"]
    timestamps= []
    values = []
    try:
        if cycletime in time_periods:
            regex = "DWHST_TEST::R*HD1:CURRENT:" + cycletime + "MS"
            pv_names = get_pv_names(search_pvs(regex, pvs_addr=pvs_addr))
            for i in pv_names:
                tries = 1
                while tries <= max_tries:
                    try:
                        ctxt.get(i)
                        timestamps.append(datetime.datetime.fromtimestamp(ctxt.get(i).timestamp))
                        values.append(str(ctxt.get(i).raw["value"]))
                        break
                    except Exception as e:
                        print(f"Got {e}. Trying again ({max_tries - tries} left).")
                    finally:
                            tries += 1  
                if tries > max_tries:
                    print("Ran out of tries. Moving onto the next value.")    
        else:
            print("Incorrect cycletime entered.")
    except Exception as e:
        print(f"Got {e}. {'This means that it could not connect to the server' if e == TimeoutError else 'Unknown cause.'}")
    final_list = [pv_names, values, timestamps]
    
    if save: save_to_csv(convert_to_df(final_list), filename)
    return convert_to_df(final_list)

def get_EPICS_Vertical_Correctors(cycletime: str, filename: str | None = None, save: bool = False, pvs_addr: str = "http://athena.isis.rl.ac.uk:9505/getPVStatus", live_addr: str = "infra.isis.rl.ac.uk:7075", max_tries: int = 3) -> DataFrame: #in MS
    """Function to get live Vertical Corrector data from EPICS. Returns a DataFrame, given the time in cycle. Will save to filename if save is True. If filename is None, will automatically set the filename.  This function is dependent on p4p, see module description. This function will not work unless p4p imported."""
    os.environ["EPICS_PVA_NAME_SERVERS"] = live_addr
    ctxt = Context("pva")
    time_periods = ["-.4", "-.2", "0", ".5", "1", "1.5", "2", "2.5", "3", "3.5", "4", "4.5", "5", "5.5", "6", "6.5", "7", "7.5", "8.1", "9", "9.8", "10.1"]
    timestamps= []
    values = []
    if filename is None: filename='EPICS_Tune_VC_'+str(cycletime)+'.dat'
    try:
        if cycletime in time_periods:
            regex = "DWVST_TEST::R*VD1:CURRENT:" + cycletime + "MS"
            pv_names = get_pv_names(search_pvs(regex, pvs_addr=pvs_addr))
            for i in pv_names:
                tries = 1
                while tries <= max_tries:
                    try:
                        timestamp = datetime.datetime.fromtimestamp(ctxt.get(i).timestamp).replace(microsecond=0)
                        values.append(str(ctxt.get(i).raw["value"]))
                        break
                    except Exception as e:
                        print(f"Got {e}. Trying again ({max_tries - tries} left).")
                    finally:
                            tries += 1  
                if tries > max_tries:
                    print("Ran out of tries. Moving onto the next value.")
        else:
            print("Incorrect cycletime entered.")
    except Exception as e:
        print(f"Got {e}. {'This means that it could not connect to the server' if e == TimeoutError else 'Unknown cause.'}")
    final_list = [pv_names, values, timestamps]
    
    if save: save_to_csv(convert_to_df(final_list), filename)
    return convert_to_df(final_list)

def get_EPICS_Tune(cycletime: str, filename: str | None = None, save: bool = False, pvs_addr: str = "http://athena.isis.rl.ac.uk:9505/getPVStatus", live_addr: str = "infra.isis.rl.ac.uk:7075", max_tries: int = 3) -> DataFrame: #in MS
    """Function to get live Tune data from EPICS. Returns a DataFrame, given the time in cycle. Will save to filename if save is True. If filename is None, will automatically set the filename.  This function is dependent on p4p, see module description. This function will not work unless p4p imported."""

    os.environ["EPICS_PVA_NAME_SERVERS"] = live_addr
    ctxt = Context("pva")
    
    if filename is None: filename='EPICS_Tune_'+str(cycletime)+'.dat'
        
    time_periods = ["-.6", "-.5", "-.4", "-.3", "-.2", "-.1", "0", ".5", "1", "1.5", "2", "2.5", "3", "3.5", "4", "4.5", "5", "5.5", "6", "7", "8", "9", "10"]
    timestamps = []
    values = []
    try:
        if cycletime in time_periods:
            regex = "DWTRIM::*_Q:AT_TIME:" + cycletime + "MS"
            pv_names = get_pv_names(search_pvs(regex, pvs_addr=pvs_addr))
            for i in pv_names:
                tries = 1
            
                while tries <= max_tries:
                    try:    
                            #timestamps.append(datetime.datetime.fromtimestamp(ctxt.get(i).timestamp))
                            timestamps.append(datetime.datetime.fromtimestamp(ctxt.get(i).timestamp).replace(microsecond=0))
                            values.append(str(ctxt.get(i).raw["value"]))
                            break
                    except Exception as e:
                        print(f"Got {e}. Trying again ({max_tries - tries} left).")
                    finally:
                            tries += 1  
                if tries > max_tries:
                    print("Ran out of tries. Moving onto the next value.")
        else:
            print("Incorrect cycle_time entered.")
    except Exception as e:
        print(f"Got {e}. {'This means that it could not connect to the server' if e == TimeoutError else 'Unknown cause.'}")
    final_list = [pv_names, values, timestamps]
    df = convert_to_df(final_list)
    df = df.transpose()
    df.columns = ['PV', 'Q_request', 'Last_change']
    df.reset_index(drop=True, inplace=True)
    if save: save_to_csv(df, filename)
    return df

def get_EPICS_HD(cycle_time: str, filename: str | None = None, save: bool = False, pvs_addr: str = "http://athena.isis.rl.ac.uk:9505/getPVStatus", live_addr: str = "infra.isis.rl.ac.uk:7075", debug_pv: bool = False, max_tries: int = 3) -> DataFrame: #in MS
    """Function to get live Horizontal Dipole data from EPICS. Returns a DataFrame, given the time in cycle. Will save to filename if save is True. If filename is None, will automatically set the filename.  This function is dependent on p4p, see module description. This function will not work unless p4p imported."""

    os.environ["EPICS_PVA_NAME_SERVERS"] = live_addr
    ctxt = Context("pva")
    
    cycle_time = str(cycle_time)
    
    if filename is None:
        filename = 'EPICS_HD_' + str(cycle_time) + '.dat'
        
    time_periods = ["-.4", "-.2", "0", ".5", "1", "1.5", "2", "2.5", "3", "3.5", "4", "4.5", "5", "5.5", "6", "6.5", "7", "7.5", "8", "8.5", "9", "10"]
    timestamps = []
    values = []
    correctors = []
    cycle_times = []
    try:
        if cycle_time in time_periods:
            regex = "DWHST_TEST::R*HD1:CURRENT:" + cycle_time + "MS"
            pv_names = get_pv_names(search_pvs(regex, pvs_addr=pvs_addr))
            if debug_pv:
                for i in pv_names:
                    print(i)
                
            for i in pv_names:
                tries = 1
                while tries <= max_tries:
                    try:
                        #print(ctxt.get(i))
                        timestamps.append(datetime.datetime.fromtimestamp(ctxt.get(i).timestamp).replace(microsecond=0))
                        values.append(round_sig(float(ctxt.get(i).raw["value"]), 4))
                        correctors.append(i.split("::")[1].split(":")[0])
                        cycle_times.append(float(i.split(":CURRENT:")[1].replace("MS", "")))  
                        break
                    except Exception as e:
                        print(f"Got {e}. Trying again ({max_tries - tries} left).")
                    finally:
                        tries += 1  
                if tries > max_tries:
                    print("Ran out of tries. Moving onto the next value.")
        else:
            print("Incorrect cycle_time entered.")
            return None
    except Exception as e:
        print(f"Got {type(e)}. {'This means that it could not connect to the server' if e == TimeoutError else 'Unknown cause.'}")
    final_list = [correctors, cycle_times, values, timestamps, pv_names]
    
    df = pd.DataFrame(final_list).transpose()
    df.columns = ['Corrector', 'Cycle_Time', 'Current', 'Last_change', 'PV']
    df.reset_index(drop=True, inplace=True)
    
    if save: save_to_csv(df, filename)
    return df

def get_EPICS_VD(cycle_time: str, filename: str | None = None, save: bool = False, pvs_addr: str = "http://athena.isis.rl.ac.uk:9505/getPVStatus", live_addr: str = "infra.isis.rl.ac.uk:7075", max_tries: int = 3) -> DataFrame: #in MS
    """Function to get live Vertical Dipole data from EPICS. Returns a DataFrame, given the time in cycle. Will save to filename if save is True. If filename is None, will automatically set the filename.  This function is dependent on p4p, see module description. This function will not work unless p4p imported."""

    os.environ["EPICS_PVA_NAME_SERVERS"] = live_addr
    ctxt = Context("pva")
    
    cycle_time = str(cycle_time)
    
    if filename is None:
        filename = 'EPICS_VD_' + str(cycle_time) + '.dat'
        
    time_periods = ["-.4", "-.2", "0", ".5", "1", "1.5", "2", "2.5", "3", "3.5", "4", "4.5", "5", "5.5", "6", "6.5", "7", "7.5", "8.1", "9", "9.8", "10.1"]
    timestamps = []
    values = []
    correctors = []
    cycle_times = []
    try:
        if cycle_time in time_periods:
            regex = "DWVST_TEST::R*VD1:CURRENT:" + cycle_time + "MS"
            pv_names = get_pv_names(search_pvs(regex, pvs_addr=pvs_addr))
            for i in pv_names:
                tries = 1
                while tries <= max_tries:
                    try:
                        timestamps.append(datetime.datetime.fromtimestamp(ctxt.get(i).timestamp).replace(microsecond=0))
                        values.append(round_sig(float(ctxt.get(i).raw["value"]), 4))
                        correctors.append(i.split("::")[1].split(":")[0])
                        cycle_times.append(float(i.split(":CURRENT:")[1].replace("MS", "")))
                        break
                    except Exception as e:
                        print(f"Got {e}. Trying again ({max_tries - tries} left).")
                    finally:
                            tries += 1  
                if tries > max_tries:
                    print("Ran out of tries. Moving onto the next value.")
        else:
            print("Incorrect cycle_time entered.")
            return None
    except Exception as e:
        print(f"Got {e}. {'This means that it could not connect to the server' if e == TimeoutError else 'Unknown cause.'}")
    final_list = [correctors, cycle_times, values, timestamps, pv_names]
    
    df = pd.DataFrame(final_list).transpose()
    df.columns = ['Corrector', 'Cycle_Time', 'Current', 'Last_change', 'PV']
    df.reset_index(drop=True, inplace=True)
    
    if save: save_to_csv(df, filename)
    return df

def get_EPICS_Q_full_cycle(filename: str | None = None, save: bool = False, pvs_addr: str = "http://athena.isis.rl.ac.uk:9505/getPVStatus", live_addr: str = "infra.isis.rl.ac.uk:7075") -> DataFrame: 
    """Function to get live Q/Tune data from EPICS. Returns a DataFrame, for the full cycle. Will save to filename if save is True. If filename is None, will automatically set the filename.  This function is dependent on p4p, see module description. This function will not work unless p4p imported."""

    time_periods = ["0", ".5", "1", "1.5", "2", "2.5", "3", "3.5", "4", "5", "5.5", "6", "7", "8", "10"]
    
    # Initialize an empty list to store the DataFrames
    df_list = []
    
    if filename is None:
        filename = 'EPICS_Q_FULL_CYCLE_.dat'
    # Iterate through each cycle time and collect the DataFrames
    for cycle_time in time_periods:
        try:
            df = get_EPICS_Tune(cycle_time, pvs_addr=pvs_addr, live_addr=live_addr) # in MS
            df_list.append(df)
        except Exception as e:
            print(f"An error occurred for cycle time {cycle_time}: {e}")
            continue
    
    # Concatenate all DataFrames into a single DataFrame
    full_cycle_df = pd.concat(df_list, ignore_index=True)
    
    # Extract cycle_time from PV and create a new column
    full_cycle_df['cycle_time'] = full_cycle_df['PV'].str.extract(r'AT_TIME:([\d\.]+)MS')
    
    # Split the DataFrame into H_Q and V_Q
    df_hq = full_cycle_df[full_cycle_df['PV'].str.contains('H_Q')].copy()
    df_vq = full_cycle_df[full_cycle_df['PV'].str.contains('V_Q')].copy()
    
    # Rename columns for clarity
    df_hq.rename(columns={'Q_request': 'Qh', 'Last_change': 'Last_change_Qh'}, inplace=True)
    df_vq.rename(columns={'Q_request': 'Qv', 'Last_change': 'Last_change_Qv'}, inplace=True)
    
    # Merge the two DataFrames on the cycle_time
    result_df = pd.merge(df_hq[['cycle_time', 'Qh', 'Last_change_Qh']], df_vq[['cycle_time', 'Qv', 'Last_change_Qv']], on='cycle_time', how='outer')
    
    # Convert Qh and Qv columns to numeric
    result_df['Qh'] = pd.to_numeric(result_df['Qh'], errors='coerce')
    result_df['Qv'] = pd.to_numeric(result_df['Qv'], errors='coerce')
    
    # Apply round_sig to Qh and Qv columns
    result_df['Qh'] = result_df['Qh'].apply(lambda x: round_sig(x, 3) if pd.notnull(x) else x)
    result_df['Qv'] = result_df['Qv'].apply(lambda x: round_sig(x, 3) if pd.notnull(x) else x)
    
    # Rearrange the columns as requested
    result_df = result_df[['cycle_time', 'Qh', 'Qv', 'Last_change_Qh', 'Last_change_Qv']]
    
    # Optionally save the combined DataFrame to a CSV file
    
    if save: result_df.to_csv(filename, index=False)
    
    return result_df

def get_EPICS_HD_full_cycle(filename: str | None = None, save: bool = False, pvs_addr: str = "http://athena.isis.rl.ac.uk:9505/getPVStatus", live_addr: str = "infra.isis.rl.ac.uk:7075") -> DataFrame: 
    """Function to get live Horisontal Dipole data from EPICS. Returns a DataFrame, for the while cycle. Will save to filename if save is True. If filename is None, will automatically set the filename.  This function is dependent on p4p, see module description. This function will not work unless p4p imported."""

    time_periods = ["0", ".5", "1", "1.5", "2", "2.5", "3", "3.5", "4", "4.5", "5", "5.5", "6", "6.5", "7", "7.5", "8", "8.5", "9", "10"]
    
    # Initialize an empty list to store the DataFrames
    df_list = []
    

    if filename is None:
        filename = 'EPICS_HD_FULL_CYCLE_.dat'
    # Iterate through each cycle time and collect the DataFrames
    for cycle_time in time_periods:
        try:
            df = get_EPICS_HD(cycle_time, pvs_addr=pvs_addr, live_addr=live_addr) # in MS
            df_list.append(df)
            #print(df)
        except Exception as e:
            print(f"An error occurred for cycle time {cycle_time}: {e}")
            continue
    
    # Concatenate all DataFrames into a single DataFrame
    full_cycle_df = pd.concat(df_list, ignore_index=True)
    
    
    # Optionally save the combined DataFrame to a CSV file
    
    if save: full_cycle_df.to_csv(filename, index=False)
    
    return full_cycle_df

def get_EPICS_VD_full_cycle(filename: str | None = None, save: bool = False, pvs_addr: str = "http://athena.isis.rl.ac.uk:9505/getPVStatus", live_addr: str = "infra.isis.rl.ac.uk:7075") -> DataFrame: 
    """Function to get live Horisontal Dipole data from EPICS. Returns a DataFrame, for the while cycle. Will save to filename if save is True. If filename is None, will automatically set the filename.  This function is dependent on p4p, see module description. This function will not work unless p4p imported."""

    time_periods = ["0", ".5", "1", "1.5", "2", "2.5", "3", "3.5", "4", "4.5", "5", "5.5", "6", "6.5", "7", "7.5", "8", "8.5", "9", "10"]
    
    # Initialize an empty list to store the DataFrames
    df_list = []
    

    if filename is None:
        filename = 'EPICS_VD_FULL_CYCLE_.dat'
    # Iterate through each cycle time and collect the DataFrames
    for cycle_time in time_periods:
        try:
            df = get_EPICS_VD(cycle_time, pvs_addr=pvs_addr, live_addr=live_addr) # in MS
            df_list.append(df)
        except Exception as e:
            print(f"An error occurred for cycle time {cycle_time}: {e}")
            continue
    
    # Concatenate all DataFrames into a single DataFrame
    full_cycle_df = pd.concat(df_list, ignore_index=True)
    
    
    # Optionally save the combined DataFrame to a CSV file
    
    if save: full_cycle_df.to_csv(filename, index=False)
    
    return full_cycle_df

def get_EPICS_Harmonics(cycletime: str, filename: str | None = None, save: bool = False, pvs_addr: str = "http://athena.isis.rl.ac.uk:9505/getPVStatus", live_addr: str = "infra.isis.rl.ac.uk:7075", max_tries: int = 3) -> DataFrame: #in MS
    """Function to get live Harmonics data from EPICS. Returns a DataFrame, given the time in cycle. Will save to filename if save is True. If filename is None, will automatically set the filename.  This function is dependent on p4p, see module description. This function will not work unless p4p imported."""

    os.environ["EPICS_PVA_NAME_SERVERS"] = live_addr
    ctxt = Context("pva")
    
    if filename is None: filename='EPICS_Harmonics_'+str(cycletime)+'.dat'
        
    time_periods = ["-.6", "-.5", "-.4", "-.3", "-.2", "-.1", "0", ".5", "1", "1.5", "2", "2.5", "3", "3.5", "4", "4.5", "5", "5.5", "6", "7", "8", "9", "10"]
    timestamps = []
    values = []
    try:
        if cycletime in time_periods:
            regex = "DWTRIM::*COS:AT_TIME:" + cycletime + "MS"
            pv_names = get_pv_names(search_pvs(regex, pvs_addr=pvs_addr))
            regex = "DWTRIM::*SIN:AT_TIME:" + cycletime + "MS"
            for n in get_pv_names(search_pvs(regex, pvs_addr=pvs_addr)): pv_names.append(n)
            for i in pv_names:
                print(i)
                tries = 1
            
                while tries <= max_tries:
                    try:    
                            #timestamps.append(datetime.datetime.fromtimestamp(ctxt.get(i).timestamp))
                            timestamps.append(datetime.datetime.fromtimestamp(ctxt.get(i).timestamp).replace(microsecond=0))
                            values.append(str(ctxt.get(i).raw["value"]))
                            break
                    except Exception as e:
                        print(f"Got {e}. Trying again ({max_tries - tries} left).")
                    finally:
                            tries += 1  
                if tries > max_tries:
                    print("Ran out of tries. Moving onto the next value.")
        else:
            print("Incorrect cycle_time entered.")
    except Exception as e:
        print(f"Got {e}. {'This means that it could not connect to the server' if e == TimeoutError else 'Unknown cause.'}")
    final_list = [pv_names, values, timestamps]
    df = convert_to_df(final_list)
    df = df.transpose()
    df.columns = ['PV', 'Harmonic', 'Last_change']
    df.reset_index(drop=True, inplace=True)
    if save: save_to_csv(df, filename)
    return df

def get_EPICS_Harmonics_full_cycle(filename: str | None = None, save: bool = False, pvs_addr: str = "http://athena.isis.rl.ac.uk:9505/getPVStatus", live_addr: str = "infra.isis.rl.ac.uk:7075", max_tries: int = 3) -> DataFrame: #in MS
    """Function to get live Harmonics data from EPICS. Returns a DataFrame, for the while cycle. Will save to filename if save is True. If filename is None, will automatically set the filename.  This function is dependent on p4p, see module description. This function will not work unless p4p imported."""

    time_periods = ["0", ".5", "1", "1.5", "2", "2.5", "3", "3.5", "4", "4.5", "5", "5.5", "6", "6.5", "7", "7.5", "8", "8.5", "9", "10"]
    
    # Initialize an empty list to store the DataFrames
    df_list = []
    

    if filename is None:
        filename = 'EPICS_Harmonics_FULL_CYCLE_.dat'
    # Iterate through each cycle time and collect the DataFrames
    for cycle_time in time_periods:
        try:
            df = get_EPICS_Harmonics(cycle_time, pvs_addr=pvs_addr, live_addr=live_addr) # in MS
            df_list.append(df)
        except Exception as e:
            print(f"An error occurred for cycle time {cycle_time}: {e}")
            continue
    
    # Concatenate all DataFrames into a single DataFrame
    full_cycle_df = pd.concat(df_list, ignore_index=True)
    
    
    # Optionally save the combined DataFrame to a CSV file
    
    if save: full_cycle_df.to_csv(filename, index=False)
    
    return full_cycle_df