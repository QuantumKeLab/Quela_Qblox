""" Go qblox_drive_AS.Configs.Manuall_QG_manage.py set your dressed state readout attenuation first """
from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager
from qblox_drive_AS.support.ExpFrames import Dressed_CavitySearching
#// test okay.


''' fill in '''
Execution:bool = True
DRandIP = {"dr":"drke","last_ip":"242"}
freq_range:dict = { #FQV1
                  "q0":[4.8085e9, 4.811e9],
                   "q1":[4.91e9, 4.913e9],
                   "q2":[5.009e9, 5.0115e9],
                   "q3":[5.11e9, 5.114e9]
                }   # np.linspace(rof+span, rof+span, freq_pts)
ro_amp = {"q0":0.07, "q1":0.02, "q2":0.08, "q3":0.04}#, "q1":1

freq_pts:int = 300
AVG:int = 200

''' Don't Touch '''
save_dir = Data_manager().build_packs_folder()
EXP = Dressed_CavitySearching(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),data_folder=save_dir)
EXP.SetParameters(freq_range,ro_amp,freq_pts,AVG,Execution)
EXP.WorkFlow()
EXP.RunAnalysis()