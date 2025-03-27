from qblox_drive_AS.support.Path_Book import find_latest_QD_pkl_for_dr
from qblox_drive_AS.support import Data_manager
from qblox_drive_AS.support.ExpFrames import QubitMonitor
#// test okay

DRandIP = {"dr":"drke","last_ip":"242"}

QM = QubitMonitor(QD_path=find_latest_QD_pkl_for_dr(DRandIP["dr"],DRandIP["last_ip"]),save_dir=Data_manager().build_packs_folder(),execution=True)


# T1 settings
QM.T1_time_range = {"q3":[0,1000e-6]} # skip if empty.
# T2 settings
QM.T2_time_range = {"q3":[0,30e-6]} # skip if empty.
QM.echo_pi_num = 0                                # int, if 0 ramsey, 1 spin echo, more for CPMG. All qubit in T2_time_range share a same number.
QM.a_little_detune_Hz = 0.2e6                     # for all qubit in T2_time_range. Set None if you do SpinEcho or CPMG
# SingleShot settings, skip if 0 shot
QM.OS_shots = 10000     
QM.OS_target_qs = ["q3"]


# T1, T2 shared settings
QM.time_sampling_func = "linspace"
QM.time_ptsORstep = 100
QM.AVG = 500


""" While looping start """
QM.StartMonitoring()
