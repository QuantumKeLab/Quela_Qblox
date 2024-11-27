from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager
from qblox_drive_AS.support.ExpFrames import DragCali
#// testing...


''' fill in '''
Execution:bool = 1
DRandIP = {"dr":"dr1","last_ip":"11"}
drag_coef_range:dict = {"q4":[-2,0]}
coef_sampling_func:str = 'linspace'
coef_ptsORstep:int|float = 50
AVG:int = 500
seq_repeat_N:int = 20

''' Don't Touch '''
save_dir = Data_manager().build_packs_folder()
EXP = DragCali(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),data_folder=save_dir)
EXP.SetParameters(drag_coef_range,coef_sampling_func,coef_ptsORstep,seq_repeat_N,avg_n=AVG,execution=Execution)
EXP.WorkFlow()
EXP.RunAnalysis()