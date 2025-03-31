from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager
from qblox_drive_AS.support.ExpFrames import Ramsey
#// test okay

''' fill in '''
Execution:bool = 1
DRandIP = {"dr":"drke","last_ip":"242"}
time_range:dict = {"q3":[0,15e-6]}
time_sampling_func:str = "linspace"
time_ptsORstep:int|float = 100
AVG:int = 500
histo_counts:int = 30

#?? Notes: While you're using one-shot method, time samples (time_pts) * shots (AVG) must less than or equal to 131000. 
use_OneShot:bool = False

''' Don't Touch '''
save_dir = Data_manager().build_packs_folder()
EXP = Ramsey(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),data_folder=save_dir)
EXP.SetParameters(time_range,time_sampling_func,time_ptsORstep,histo_counts,AVG,Execution,OSmode=use_OneShot)
EXP.WorkFlow()
EXP.RunAnalysis()
