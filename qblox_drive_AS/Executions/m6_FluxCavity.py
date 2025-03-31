from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager
from qblox_drive_AS.support.ExpFrames import FluxCavity
#// test okay.

''' fill in '''
Execution:bool = True
DRandIP = {"dr":"drke","last_ip":"242"}
freq_span_range:dict = {"q3":[-3e6,+3e6]}#"q0":[-3e6,+3e6], "q1":[-2e6,+2e6],"q2":[-2.5e6,+2.5e6], "q3":[-3e6,+3e6]}    # np.linspace(rof+span, rof+span, freq_pts)
flux_range:list = [-1.75, 1.75, 40]                                 # flux [from, end, pts/step]
flux_sampling_func:str = 'linspace'                          # 'linspace'/ 'logspace'/ 'arange

freq_pts:int = 40
AVG:int = 100

''' Don't Touch '''
save_dir = Data_manager().build_packs_folder()
EXP = FluxCavity(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),data_folder=save_dir)
EXP.SetParameters(freq_span_range,flux_range,flux_sampling_func,freq_pts,AVG,Execution)
EXP.WorkFlow()
EXP.RunAnalysis()
