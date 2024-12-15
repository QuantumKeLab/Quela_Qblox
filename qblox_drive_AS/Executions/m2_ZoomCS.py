from qblox_drive_AS.support.ExpFrames import Zoom_CavitySearching
from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager
#// Test okay.

''' fill in '''
Execution:bool = True
DRandIP = {"dr":"drke","last_ip":"242"}
freq_range:dict = {
                  #FQV1
                  # "q0":[4.505e9, 4.525e9],
                  #  "q1":[4.59e9, 4.607e9],
                  #  "q2":[4.690e9, 4.72e9],
                  #  "q3":[4.795e9, 4.804e9],
 
                
                #5Q array
                "q0":[5.8450e9, 5.8550e9],
                  "q1":[6.1e9, 6.13e9],
                  "q2":[5.92e9, 5.93e9],
                  "q3":[5.72e9, 5.74e9],
                  # "q0":[5.80e9, 5.9e9],
                  # "q1":[6.1e9, 6.2e9],
                   }
freq_pts:int = 100
AVG:int = 100

''' Don't Touch '''
save_dir = Data_manager().build_packs_folder()
EXP = Zoom_CavitySearching(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),data_folder=save_dir)
EXP.SetParameters(freq_range,freq_pts,AVG,Execution)
EXP.WorkFlow()
EXP.RunAnalysis()

