states = ["zero"]
shotnums = [500, 1000, 2000, 5000, 10000]
# fre = ["09995","09996","09997","09998","09999","10000","10001","10002","10003","10004","10005"]
fre = ["50","100","200","400","-50","-100","-200","-400"]
amp = ["0000","0125","0250","0375","0500","0625","0750","0875","1000"]
from qcodes import Instrument
# for l in range(1):
#     for m in amp:
#         from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
#         from qblox_drive_AS.support import Data_manager
#         from qblox_drive_AS.support.ExpFrames import XTomography
#         #// test okay.

#         ''' fill in '''
#         Execution: bool = 1
#         DRandIP = {"dr": "dr1", "last_ip": "11"}
#         target_q: list = ["q2"]
#         control_q: str = "q2"
#         noise_q: str = "q0"
#         un_trained_pulse: bool = 1
#         maxmum_gate = 101
#         avg_n: int = 5000
#         init_state: str = "zero" # zero plus plus_i
#         ''' Don't Touch '''
#         gate_counts = range(maxmum_gate)
#         save_dir = Data_manager().build_packs_folder()
#         EXP = XTomography(QD_path=f"/Users/kelab_01/Documents/GitHub/Quela_Qblox/qblox_drive_AS/QD_backup/20260405/DR1#11_SumInfo{m}amp.pkl", data_folder=save_dir)
#         EXP.SetParameters(target_q, init_state, avg_n, gate_counts, Execution, control_q, noise_q)
#         EXP.WorkFlow()
#         EXP.RunAnalysis(m)
#         Instrument.close_all()

# for l in range(5):
#     from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
#     from qblox_drive_AS.support import Data_manager
#     from qblox_drive_AS.support.ExpFrames import XTomography
#     #// test okay.
#     realzero = "none"
#     ''' fill in '''
#     Execution: bool = 1
#     DRandIP = {"dr": "dr1", "last_ip": "11"}
#     target_q: list = ["q2"]
#     control_q: str = "q2"
#     noise_q: str = None
#     un_trained_pulse: bool = 1
#     maxmum_gate = 101
#     avg_n: int = 5000
#     init_state: str = "zero" # zero plus plus_i
#     ''' Don't Touch '''
#     gate_counts = range(maxmum_gate)
#     save_dir = Data_manager().build_packs_folder()
#     EXP = XTomography(QD_path=f"/Users/kelab_01/Documents/GitHub/Quela_Qblox/qblox_drive_AS/QD_backup/20260405/DR1#11_SumInfo0000amp.pkl", data_folder=save_dir)
#     EXP.SetParameters(target_q, init_state, avg_n, gate_counts, Execution, control_q, noise_q)
#     EXP.WorkFlow()
#     EXP.RunAnalysis(realzero)
#     Instrument.close_all()


for l in range(5):
    from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
    from qblox_drive_AS.support import Data_manager
    from qblox_drive_AS.support.ExpFrames import XTomography
    #// test okay.
    realzero = "100ns"
    ''' fill in '''
    Execution: bool = 1
    DRandIP = {"dr": "dr1", "last_ip": "11"}
    target_q: list = ["q2"]
    control_q: str = "q2"
    noise_q: str = "q0"
    un_trained_pulse: bool = 1
    maxmum_gate = 101
    avg_n: int = 5000
    init_state: str = "zero" # zero plus plus_i
    ''' Don't Touch '''
    gate_counts = range(maxmum_gate)
    save_dir = Data_manager().build_packs_folder()
    EXP = XTomography(QD_path=f"/Users/kelab_01/Documents/GitHub/Quela_Qblox/qblox_drive_AS/QD_backup/20260405/DR1#11_SumInfo100ns.pkl", data_folder=save_dir)
    EXP.SetParameters(target_q, init_state, avg_n, gate_counts, Execution, control_q, noise_q)
    EXP.WorkFlow()
    EXP.RunAnalysis(realzero)
    Instrument.close_all()


# for k in shotnums:
#     for j in range(min(8, int(30000/k))):
#         for i in states:
#             from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
#             from qblox_drive_AS.support import Data_manager
#             from qblox_drive_AS.support.ExpFrames import XTomography
#             #// test okay.

#             ''' fill in '''
#             Execution: bool = 1
#             DRandIP = {"dr": "dr1", "last_ip": "11"}
#             target_q: list = ["q3"]
#             control_q: str = "q3"
#             noise_q: str = None
#             un_trained_pulse: bool = True
#             maxmum_gate = 101
#             avg_n: int = k
#             init_state: str = i # zero plus plus_i
#             ''' Don't Touch '''
#             gate_counts = range(maxmum_gate)
#             save_dir = Data_manager().build_packs_folder()
#             EXP = XTomography(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"], DRandIP["last_ip"]), data_folder=save_dir)
#             EXP.SetParameters(target_q, init_state, avg_n, gate_counts, Execution, control_q, noise_q)
#             EXP.WorkFlow()
#             EXP.RunAnalysis(i,k)
#             Instrument.close_all()