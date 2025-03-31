""" Go qblox_drive_AS.Configs.Manuall_QG_manage.py set your dressed state readout attenuation first """
from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager
from qblox_drive_AS.support.ExpFrames import Dressed_CavitySearching
#// test okay.


''' fill in '''
Execution:bool = True
DRandIP = {"dr":"drke","last_ip":"242"}
freq_range:dict = {"q0":[4.8045e9, 4.8145e9],
                   "q1":[4.905e9, 4.915e9],
                   "q2":[5.005e9, 5.015e9],
                   "q3":[5.105e9, 5.120e9],}    # np.linspace(rof+span, rof+span, freq_pts)
ro_amp = {"q0":0.1, "q1":0.05,"q2":0.12, "q3":0.05}

freq_pts:int = 500
AVG:int = 500

''' Don't Touch '''
save_dir = Data_manager().build_packs_folder()
EXP = Dressed_CavitySearching(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),data_folder=save_dir)
EXP.SetParameters(freq_range,ro_amp,freq_pts,AVG,Execution)
EXP.WorkFlow()
EXP.RunAnalysis()