""" Go qblox_drive_AS.Configs.Manuall_QG_manage.py set your dressed state readout attenuation first """
from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager
from qblox_drive_AS.support.ExpFrames import Dressed_CavitySearching
#// test okay.


''' fill in '''
Execution:bool = True
DRandIP = {"dr":"drke","last_ip":"242"}
freq_range:dict = {
                #FQV1
                # "q0":[4.505e9, 4.525e9], "q1":[4.5995e9, 4.6015e9],
                # "q2":[4.690e9, 4.72e9], "q3":[4.795e9, 4.804e9]
                   
                   #5Q array
                   "q0":[5.8450e9, 5.8550e9], "q1":[6.115e9, 6.123e9],
                   "q2":[5.90e9, 5.94e9],"q3":[5.72e9, 5.74e9]}    # np.linspace(rof+span, rof+span, freq_pts)
ro_amp = {"q0":0.2, "q1":0.2,"q2":0.1, "q3":0.1}#, "q3":0.05}#"q1":0.1, 
freq_pts:int = 100
AVG:int = 100

''' Don't Touch '''
save_dir = Data_manager().build_packs_folder()
EXP = Dressed_CavitySearching(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),data_folder=save_dir)
EXP.SetParameters(freq_range,ro_amp,freq_pts,AVG,Execution)
EXP.WorkFlow()
EXP.RunAnalysis()