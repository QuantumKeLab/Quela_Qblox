import os, sys
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from Modularize.support import QDmanager
from Modularize.support.Path_Book import find_latest_QD_pkl_for_dr


def fillin_PDans(QD_path:str,ans:dict):
    """
    Fill in the power dep answer to the quantum_device.\n
    format:\n
    `ans = {"q0":{"dressF_Hz":,"dressP":,"bareF_Hz":},...}`
    """
    QDagent = QDmanager(QD_path)
    QDagent.QD_loader()
    for q in ans:
        qubit = QDagent.quantum_device.get_element(q)
        if ans[q]["dressP"] != "": qubit.measure.pulse_amp(ans[q]["dressP"]) 
        if ans[q]["dressF_Hz"] != "": qubit.clock_freqs.readout(ans[q]["dressF_Hz"])
        if ans[q]["bareF_Hz"] != "": QDagent.Notewriter.save_bareFreq_for(target_q=q,bare_freq=ans[q]["bareF_Hz"])
        if ans[q]["ro_atte"] != "": QDagent.Notewriter.save_DigiAtte_For(atte_dB=ans[q]["ro_atte"],target_q=q,mode='ro')

    QDagent.refresh_log("PD answer stored!")
    QDagent.QD_keeper()


if __name__ == "__main__":

    """ Fill in """
    DRandIP = {"dr":"drke","last_ip":"116"}
    PDans = {
        "q0":{"dressF_Hz":6.11077e9,"dressP":0.7,"bareF_Hz":6.100552e9,"ro_atte":50}, #q4
        # "q1":{"dressF_Hz":5.934e9,"dressP":0.3,"bareF_Hz":5.9206e9,"ro_atte":40},   #q3
        # "q1":{"dressF_Hz":6.0775e9,"dressP":0.6,"bareF_Hz":6.070e9,"ro_atte":50},   #q2
        # "q1":{"dressP":0.2,"ro_atte":40},   #q2


    }
    

    """ Storing """
    QD_path = find_latest_QD_pkl_for_dr(which_dr=DRandIP["dr"],ip_label=DRandIP["last_ip"])
    fillin_PDans(QD_path, PDans)

    
    
