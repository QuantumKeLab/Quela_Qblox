import os, sys 
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', ".."))
from qblox_drive_AS.support.Pulse_schedule_library import mix_T1_sche, T1_sche, set_LO_frequency, pulse_preview, IQ_data_dis, dataset_to_array, T1_fit_analysis, Fit_analysis_plot
from numpy import median, std, array
from qblox_drive_AS.support import QDmanager, Data_manager
from xarray import open_dataset

def hito_replot(folder_path,QD_agent,qubit):
    t1_us_rec = []
    files = [os.path.join(folder_path,name) for name in os.listdir(folder_path) if (os.path.isfile(os.path.join(folder_path,name)) and name.split(".")[-1]=="nc")]
    ref_IQ = QD_agent.refIQ[qubit]
    for file in files:
        ds = open_dataset(file)
        I,Q= dataset_to_array(dataset=ds,dims=1)
        data = IQ_data_dis(I,Q,ref_I=ref_IQ[0],ref_Q=ref_IQ[-1])
        samples = array(ds.x0)
        
        data_fit= T1_fit_analysis(data=data,freeDu=samples,T1_guess=20e-6)
        t1_us_rec.append(data_fit.attrs['T1_fit']*1e6)

    Data_manager().save_histo_pic(QD_agent,{str(qubit):t1_us_rec},qubit,mode="t1")

def meas_result_plot(file_path,QD_agent,qubit):
    ref_IQ = QD_agent.refIQ[qubit]
    ds = open_dataset(file_path)
    I,Q= dataset_to_array(dataset=ds,dims=1)
    data = IQ_data_dis(I,Q,ref_I=ref_IQ[0],ref_Q=ref_IQ[-1])
    samples = array(ds.x0)
    
    data_fit= T1_fit_analysis(data=data,freeDu=samples,T1_guess=70e-6)
    Fit_analysis_plot(data_fit,P_rescale=False,Dis=None)

# if __name__ == '__main__':
#     qubit = 'q3'
#     path = r'Modularize\Meas_raw\2024_9_9\DR4q0_T1(31)_H19M3S22.nc'
#     QD_path = r'Modularize\QD_backup\2024_9_9\DR4#81_SumInfo.pkl'
#     QD_agent = QDmanager(QD_path)
#     QD_agent.QD_loader()
#     meas_result_plot(path,QD_agent,qubit)


if __name__ == '__main__':
    qubit = 'q3'
    folder_path = r'Modularize\Meas_raw\T1'
    QD_path = r'Modularize\QD_backup\2024_9_11\DR4#81_SumInfo.pkl'
    QD_agent = QDmanager(QD_path)
    QD_agent.QD_loader()
    hito_replot(folder_path,QD_agent,qubit)