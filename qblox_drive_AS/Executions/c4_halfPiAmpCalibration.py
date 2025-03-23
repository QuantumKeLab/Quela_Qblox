from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager
from qblox_drive_AS.support.ExpFrames import hPiAcali
#// 0.9.2 okay


''' fill in '''
Execution:bool = 1
DRandIP = {"dr":"dr4","last_ip":"81"}
ro_power_coef_range:dict = {"q0":[0.9,1.1],"q1":[0.9,1.1],"q2":[0.9,1.1]}    # half pi pulse coef, it should be around 1
coef_sampling_func:str = 'linspace'
half_pi_quadruple_num:list = [3,5]
coef_ptsORstep:int|float = 50
AVG:int = 300



''' Don't Touch '''
save_dir = Data_manager().build_packs_folder()
EXP = hPiAcali(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),data_folder=save_dir)
EXP.SetParameters(ro_power_coef_range,coef_sampling_func,coef_ptsORstep,half_pi_quadruple_num,avg_n=AVG,execution=Execution)
EXP.WorkFlow()
EXP.RunAnalysis()