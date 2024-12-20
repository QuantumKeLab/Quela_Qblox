from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager
from qblox_drive_AS.support.ExpFrames import StarkShift


''' fill in '''
Execution:bool = True
DRandIP = {"dr":"drke","last_ip":"242"}
target_qs:list =['q0']
ro_p_max:float=0.12#the output figure will show (ro_p_max)**2
p_pts:int=4
freq_span_Hz:float = 500e6
AVG:int = 100

''' Don't Touch '''
ro_amp_range:list = [0, ro_p_max,p_pts]    
save_dir = Data_manager().build_packs_folder()
EXP = StarkShift(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),data_folder=save_dir)
EXP.SetParameters(target_qs, freq_span_range,ro_amp_range,ro_amp_sampling_func,freq_pts,AVG,Execution)
EXP.WorkFlow()
EXP.RunAnalysis()