from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager
from qblox_drive_AS.support.ExpFrames import ROFcali
#// test okay


''' fill in '''
Execution:bool = True
DRandIP = {"dr":"drke","last_ip":"242"}
freq_span_range:dict = {"q3":[-1e6,+1e6]}
freq_pts = 100
AVG:int = 300



''' Don't Touch '''
save_dir = Data_manager().build_packs_folder()
EXP = ROFcali(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),data_folder=save_dir)
EXP.SetParameters(freq_span_range,freq_pts,avg_n=AVG,execution=Execution)
EXP.WorkFlow()
EXP.RunAnalysis()