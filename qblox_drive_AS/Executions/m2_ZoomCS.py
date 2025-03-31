from qblox_drive_AS.support.ExpFrames import Zoom_CavitySearching
from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager
#// Test okay.

''' fill in '''
Execution:bool = True
DRandIP = {"dr":"drke","last_ip":"242"}
freq_range:dict = {"q0":[4.805e9, 4.814e9],
                   "q1":[4.9075e9, 4.914e9],
                   "q2":[5.005e9, 5.014e9],
                   "q3":[5.106e9, 5.115e9],
                #    "q4":[6.315e9, 6.325e9],
                #    "q5":[6.414e9, 6.424e9],
                }
freq_pts:int = 500
AVG:int = 500

''' Don't Touch '''
save_dir = Data_manager().build_packs_folder()
EXP = Zoom_CavitySearching(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),data_folder=save_dir)
EXP.SetParameters(freq_range,freq_pts,AVG,Execution)
EXP.WorkFlow()
EXP.RunAnalysis()

