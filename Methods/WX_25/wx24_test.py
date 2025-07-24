import os
os.environ["ACCESS_LIVE_DATA"] = "True"
os.environ["SUPPRESS_P4P_WARNING"] = "True"
import wx24
import datetime

### The following works
data = wx24.get_archive_trim_quads()
data.to_csv("get_archive_trim_quads.dat")
print(data)

data = wx24.get_archive_correctors(1, "H", datetime.datetime(2023, 5, 5, 5, 5, 5, 5), datetime.datetime(2025, 5, 5, 5, 0, 0, 0))
data.to_csv("get_archive_correctors.dat")
print(data)

data = wx24.get_EPICS_Harmonics("1")
data.to_csv("get_EPICS_Harmonics.dat")
print(data)

data = wx24.get_EPICS_HD("1")
data.to_csv("get_EPICS_HD.dat")
print(data)

data = wx24.get_EPICS_VD("1")
data.to_csv("get_EPICS_VD.dat")
print(data)

data = wx24.get_EPICS_HD_full_cycle()
data.to_csv("get_EPICS_HD_full_cycle.dat")
print(data)

data = wx24.get_EPICS_VD_full_cycle()
data.to_csv("get_EPICS_VD_full_cycle.dat")
print(data)

data = wx24.get_EPICS_Q_full_cycle()
data.to_csv("get_EPICS_Q_full_cycle.dat")
print(data)

data = wx24.get_EPICS_Horizontal_Correctors("1")
data.to_csv("get_EPICS_Horizontal_Correctors.dat")
print(data)

data = wx24.get_EPICS_Vertical_Correctors("1")
data.to_csv("get_EPICS_Vertical_Correctors.dat")
print(data)

data = wx24.get_EPICS_Harmonics_full_cycle()
data.to_csv("get_EPICS_Harmonics_full_cycle")
print(data)