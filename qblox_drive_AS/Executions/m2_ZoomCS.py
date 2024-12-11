from qblox_drive_AS.support.ExpFrames import Zoom_CavitySearching
from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager
#// Test okay.

''' fill in '''
Execution:bool = True
DRandIP = {"dr":"drke","last_ip":"242"}
freq_range:dict = {"q0":[4.50e9, 4.52e9],
                   "q1":[4.595e9, 4.615e9],
                   "q2":[4.68e9, 4.72e9],
                   "q3":[4.785e9, 4.805e9],
 
                
                #5Q array
                  # "q0":[5.8450e9, 5.8550e9],
                  # "q1":[6.1e9, 6.13e9],
                   }
freq_pts:int = 100
AVG:int = 100

''' Don't Touch '''
save_dir = Data_manager().build_packs_folder()
EXP = Zoom_CavitySearching(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),data_folder=save_dir)
EXP.SetParameters(freq_range,freq_pts,AVG,Execution)
EXP.WorkFlow()
EXP.RunAnalysis()

