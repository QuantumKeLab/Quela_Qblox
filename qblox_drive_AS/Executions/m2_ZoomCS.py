from qblox_drive_AS.support.ExpFrames import Zoom_CavitySearching
from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager
#// Test okay.

''' fill in '''
Execution:bool = True
DRandIP = {"dr":"drke","last_ip":"242"}
freq_range:dict = { #FQV1
                  "q0":[4.8080e9, 4.8105e9],
                   "q1":[4.9075e9, 4.913e9],
                   "q2":[5.0090e9, 5.011e9],
                   "q3":[5.11e9, 5.115e9]
                }
freq_pts:int = 100
AVG:int = 100

''' Don't Touch '''
save_dir = Data_manager().build_packs_folder()
EXP = Zoom_CavitySearching(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),data_folder=save_dir)
EXP.SetParameters(freq_range,freq_pts,AVG,Execution)
EXP.WorkFlow()
EXP.RunAnalysis()

