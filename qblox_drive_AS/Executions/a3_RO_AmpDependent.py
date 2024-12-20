from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager
from qblox_drive_AS.support.ExpFrames import ROAdepOS
#// test okay.

''' fill in '''
Execution:bool = True
DRandIP = {"dr":"drke","last_ip":"242"}
target_qs:list = ["q0"]#,"q1"]
roAmp_factor=[1]
#roAmp_factor=[1.3,1.4,1.5,1.6]
shots:int = 50000
histo_counts:int = 1000
''' Don't Touch '''
save_dir = Data_manager().build_packs_folder()
EXP = ROAdepOS(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),data_folder=save_dir)
EXP.SetParameters(target_qs,roAmp_factor,histo_counts,shots,Execution)
EXP.WorkFlow()
EXP.RunAnalysis(histo_ana=True if histo_counts > 1 else False)