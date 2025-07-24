import wx24
import datetime

wx24.get_EPICS_Harmonics_full_cycle().to_csv("get_EPICS_Harmonics_full_cycle.dat")
print(wx24.get_EPICS_Harmonics_full_cycle())