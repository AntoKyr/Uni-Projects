import numpy as np
from matplotlib import pyplot as plt

# Data Collection
with np.load('SimData.npz', allow_pickle=True) as SimData:
    NoR = SimData['NoR']
    phis = SimData['phis']
    Timings = SimData['Timings']
    P_data = SimData['P_data']
    T_data = SimData['T_data']
    HT_data = SimData['HT_data']
    HTC_data = SimData['HTC_data']
    Comb_data = SimData['Comb_data']
    mexh_data = SimData['mexh_data']
    mair_data = SimData['mair_data']
    mf_data = SimData['mf_data']
    int_mflow_data = SimData['int_mflow_data']
    exh_mflow_data = SimData['exh_mflow_data']
    exhint_data = SimData['exhint_data']
    V_data = SimData['V_data']
    AV_data = SimData['AV_data']
    Cf_data = SimData['Cf_data']
    ev = SimData['ev']
    exh_f = SimData['exh_f']
    imep_gross = SimData['imep_gross']
    imep_net = SimData['imep_net']
    pmep = SimData['pmep']
    effi = SimData['effi']
    dWidt = SimData['dWidt']
    comb_dur_data = SimData['comb_dur_data']
    comb_m_data = SimData['comb_m_data']
    Energy_Fracts = SimData['Energy_Fracts']

int_timings = Timings[0]
exh_timings = Timings[1]
comb_timings = Timings[2]
inj_timings = Timings[3]

int_phi_o = int_timings[0]
int_phi_c = int_timings[1]
exh_phi_o = exh_timings[0]
exh_phi_c = exh_timings[1]
comb_start = comb_timings[0]
comb_end = comb_timings[1]
inj_start = inj_timings[0]
inj_end = inj_timings[1]

def Plot_timings(lineys, Displacement=False):
    if Displacement:
        plt.plot(np.array([int_phi_o, int_phi_o]), lineys, 'go--', label='Intake')
        plt.plot(np.array([int_phi_c, int_phi_c])-720, lineys, 'go--')
        plt.plot(np.array([exh_phi_o, exh_phi_o]), lineys, 'y*--', label='Exhaust') 
        plt.plot(np.array([exh_phi_c, exh_phi_c]), lineys, 'y*--')
        plt.plot(np.array([inj_start, inj_start])-720, lineys, 'c^--', label='Injection') 
        plt.plot(np.array([inj_end, inj_end])-720, lineys, 'c^--')
        plt.plot(np.array([comb_start, comb_start])-720, lineys, 'r+--', label='Combustion') 
        plt.plot(np.array([comb_end, comb_end]), lineys, 'r+--') 
    else:
        plt.plot(np.array([int_phi_o, int_phi_o]), lineys, 'go--', label='Intake')
        plt.plot(np.array([int_phi_c, int_phi_c]), lineys, 'go--')
        plt.plot(np.array([exh_phi_o, exh_phi_o]), lineys, 'y*--', label='Exhaust') 
        plt.plot(np.array([exh_phi_c, exh_phi_c]), lineys, 'y*--')
        plt.plot(np.array([inj_start, inj_start]), lineys, 'c^--', label='Injection') 
        plt.plot(np.array([inj_end, inj_end]), lineys, 'c^--')
        plt.plot(np.array([comb_start, comb_start]), lineys, 'r+--', label='Combustion') 
        plt.plot(np.array([comb_end, comb_end]), lineys, 'r+--') 
    


# Moving some stuff around
phis2 = np.concatenate((np.array(phis[np.argwhere(phis>360)]).flatten() - 720, np.array(phis[np.argwhere(phis<360)]).flatten()))
P_plot = P_data/10**5
P_plot = np.concatenate((np.array(P_plot[np.argwhere(phis>360)]).flatten(), np.array(P_plot[np.argwhere(phis<360)]).flatten()))
HTC_plot = HTC_data
HTC_plot = np.concatenate((np.array(HTC_plot[np.argwhere(phis>360)]).flatten(), np.array(HTC_plot[np.argwhere(phis<360)]).flatten()))
T_plot = T_data
T_plot = np.concatenate((np.array(T_plot[np.argwhere(phis>360)]).flatten(), np.array(T_plot[np.argwhere(phis<360)]).flatten()))
Comb_plot = Comb_data/10**3
Comb_plot = np.concatenate((np.array(Comb_plot[np.argwhere(phis>360)]).flatten(), np.array(Comb_plot[np.argwhere(phis<360)]).flatten()))

# Pressure Diagram
plt.figure()
plt.plot(phis2, P_plot)

lineys = np.array([max(P_data), min(P_data)])/10**5
Plot_timings(lineys=lineys, Displacement=True)
plt.grid()
plt.title('Pressure Diagram')
plt.xlabel('angle [degrees]')
plt.ylabel('Pressure [bar]')
plt.legend(bbox_to_anchor=(1.1, 1.1))

# Temperature Diagram
plt.figure()
plt.plot(phis2, T_plot)

lineys = np.array([max(T_data), min(T_data)])
Plot_timings(lineys=lineys, Displacement=True)
plt.grid()
plt.title('Temperature Diagram')
plt.xlabel('angle [degrees]')
plt.ylabel('Temperature [K]')
plt.legend(bbox_to_anchor=(1.1, 1.1))

# Heat Transfer Coef Diagram
plt.figure()
plt.plot(phis2, HTC_plot)

lineys = np.array([max(HTC_data), min(HTC_data)])
Plot_timings(lineys=lineys, Displacement=True)
plt.grid()
plt.title('Heat Transfer Coef Diagram')
plt.xlabel('angle [degrees]')
plt.ylabel('Heat Transfer [W/(m**2 K)]')
plt.legend(bbox_to_anchor=(1.1, 1.1))


# Heating Power Diagram
plt.figure()
plt.plot(phis2, Comb_plot)

lineys = np.array([max(Comb_data), min(Comb_data)])/10**3
Plot_timings(lineys=lineys, Displacement=True)

plt.grid()
plt.title('Heating Power Diagram')
plt.xlabel('angle [degrees]')
plt.ylabel('Heating Power [kW]')
plt.legend(bbox_to_anchor=(1.1, 1.1))


# Mass Diagram
plt.figure()
plt.plot(phis, mexh_data*10**3, label='Exhaust Gas')
plt.plot(phis, mair_data*10**3, label='Air Mass')
plt.plot(phis, mf_data*10**3, label ='Fuel Mass')

lineys = np.array([max([max(mexh_data), max(mf_data), max(mair_data)]), min([min(mexh_data), min(mf_data), min(mair_data)])])*10**3
Plot_timings(lineys=lineys)

plt.grid()
plt.title('Mass Diagram')
plt.xlabel('angle [degrees]')
plt.ylabel('Mass [g]')
plt.legend(bbox_to_anchor=(1.1, 1.1))


# Mass Flow Diagram
plt.figure()
plt.plot(phis, int_mflow_data, label='Intake Mass Flow')
plt.plot(phis, exh_mflow_data, label='Exhaust Mass Flow')

lineys = np.array([max([max(int_mflow_data), max(exh_mflow_data)]), min([min(int_mflow_data), min(exh_mflow_data)])])
Plot_timings(lineys=lineys)

plt.grid()
plt.title('Mass Flow Diagram')
plt.xlabel('angle [degrees]')
plt.ylabel('Mass Flow Rate [kg/s]')
plt.legend(bbox_to_anchor=(1.125, 1.125))


# Exhaust Gas Settled in the Intake
plt.figure()
plt.plot(phis, exhint_data*10**3)

lineys = np.array([max(exhint_data), min(exhint_data)])*10**3
Plot_timings(lineys=lineys)

plt.grid()
plt.title('Exhaust Gas in the Intake Diagram')
plt.xlabel('angle [degrees]')
plt.ylabel('Mass [g]')
plt.legend(bbox_to_anchor=(1.1, 1.1))


# Burnt Mass Fraction Diagram
plt.figure()
plt.plot(phis, mexh_data/(mexh_data+mair_data+mf_data))

lineys = np.array([max(mexh_data/(mexh_data+mair_data+mf_data)), min(mexh_data/(mexh_data+mair_data+mf_data))])
Plot_timings(lineys=lineys)

plt.grid()
plt.title('Burnt Mass Fraction Diagram')
plt.xlabel('angle [degrees]')
plt.ylabel('f')
plt.legend(bbox_to_anchor=(1.1, 1.1))


# # Volume Diagram
# plt.figure()
# plt.plot(phis, V_data*10**3)

# lineys = np.array([max(V_data), min(V_data)])*10**3
# plt.plot(np.array([int_phi_o, int_phi_o]), lineys, 'g--', label='Intake')
# plt.plot(np.array([int_phi_c, int_phi_c]), lineys, 'g--')
# plt.plot(np.array([exh_phi_o, exh_phi_o]), lineys, 'y--', label='Exhaust') 
# plt.plot(np.array([exh_phi_c, exh_phi_c]), lineys, 'y--')
# plt.plot(np.array([inj_start, inj_start]), lineys, 'c--', label='Injection') 
# plt.plot(np.array([inj_end, inj_end]), lineys, 'c--')
# plt.plot(np.array([comb_start, comb_start]), lineys, 'r--', label='Combustion') 
# plt.plot(np.array([comb_end, comb_end]), lineys, 'r--') 

# plt.grid()
# plt.title('Volume Diagram')
# plt.xlabel('angle [degrees]')
# plt.ylabel('Volume [L]')
# plt.legend(bbox_to_anchor=(0., 1))


# # Valve Area Diagram
# plt.figure()
# plt.plot(phis, AV_data)

# lineys = np.array([max(AV_data), min(AV_data)])
# plt.plot(np.array([int_phi_o, int_phi_o]), lineys, 'g--', label='Intake')
# plt.plot(np.array([int_phi_c, int_phi_c]), lineys, 'g--')
# plt.plot(np.array([exh_phi_o, exh_phi_o]), lineys, 'y--', label='Exhaust') 
# plt.plot(np.array([exh_phi_c, exh_phi_c]), lineys, 'y--')
# plt.plot(np.array([inj_start, inj_start]), lineys, 'c--', label='Injection') 
# plt.plot(np.array([inj_end, inj_end]), lineys, 'c--')
# plt.plot(np.array([comb_start, comb_start]), lineys, 'r--', label='Combustion') 
# plt.plot(np.array([comb_end, comb_end]), lineys, 'r--') 

# plt.grid()
# plt.title('Valve Area Diagram')
# plt.xlabel('angle [degrees]')
# plt.ylabel('Area [m**2]')
# plt.legend(bbox_to_anchor=(0., 1))


# # Cf Diagram
# plt.figure()
# plt.plot(phis, Cf_data)

# lineys = np.array([max(Cf_data), min(Cf_data)])
# plt.plot(np.array([int_phi_o, int_phi_o]), lineys, 'g--', label='Intake')
# plt.plot(np.array([int_phi_c, int_phi_c]), lineys, 'g--')
# plt.plot(np.array([exh_phi_o, exh_phi_o]), lineys, 'y--', label='Exhaust') 
# plt.plot(np.array([exh_phi_c, exh_phi_c]), lineys, 'y--')
# plt.plot(np.array([inj_start, inj_start]), lineys, 'c--', label='Injection') 
# plt.plot(np.array([inj_end, inj_end]), lineys, 'c--')
# plt.plot(np.array([comb_start, comb_start]), lineys, 'r--', label='Combustion') 
# plt.plot(np.array([comb_end, comb_end]), lineys, 'r--') 

# plt.grid()
# plt.title('Cf Diagram')
# plt.xlabel('angle [degrees]')
# plt.ylabel('Cf')
# plt.legend(bbox_to_anchor=(0., 1))


# P-V Diagram
plt.figure()
plt.plot(V_data, P_data/10**5)
plt.grid()
plt.title('P-V Diagram')
plt.xlabel('V [l]')
plt.ylabel('Pressure [bar]')

# P-V Logarithmic
plt.figure()
plt.plot(V_data, P_data/10**5)
plt.grid()
plt.title('P-V Diagram')
plt.xlabel('V [l]')
plt.ylabel('Pressure [bar]')
plt.yscale('log')
plt.xscale('log')

# Convergence Diagrams

# Indicated Power 
plt.figure()
plt.plot(np.arange(0, NoR), dWidt/10**3)
plt.grid()
plt.title('Indicated Power Convergence')
plt.xlabel('No of Cycles Ran')
plt.ylabel('Power [kW]')

# Filling Efficiency
plt.figure()
plt.plot(np.arange(0, NoR), ev)
plt.grid()
plt.title('Filling Efficiency Convergence')
plt.xlabel('No of Cycles Ran')
plt.ylabel('ev')

# imep net
plt.figure()
plt.plot(np.arange(0, NoR), imep_net/10**5)
plt.grid()
plt.title('Net imep Convergence')
plt.xlabel('No of Cycles Ran')
plt.ylabel('imep [bar]')


# Indicated Efficiency
plt.figure()
plt.plot(np.arange(0, NoR), effi)
plt.grid()
plt.title('Indicated Efficiency Convergence')
plt.xlabel('No of Cycles Ran')
plt.ylabel('η')


# Combustion duration
plt.figure()
plt.plot(np.arange(0, NoR), comb_dur_data)
plt.grid()
plt.title('Combustion Duration Convergence')
plt.xlabel('No of Cycles Ran')
plt.ylabel('angle [degrees]')


# Combustion form parameter
plt.figure()
plt.plot(np.arange(0, NoR), comb_m_data)
plt.grid()
plt.title('Combustion Form Parameter Convergence')
plt.xlabel('No of Cycles Ran')
plt.ylabel('m')


# Energy Balance Pie Chart
# Qin = Energy_Fracts[-1,0]/10**3
# Qw = Energy_Fracts[-1,1]/10**3
# Wnet = Energy_Fracts[-1,2]/10**3
# Hdiff = Energy_Fracts[-1,3]/10**3
# Qdiff = Qin - Qw - Wnet - Hdiff
# Fractions = np.array([Wnet, Qw, Hdiff, Qdiff])
# print(Qin)
# print(Qw)
# print(Wnet)
# print
# plt.figure()
# plt.pie(Fractions, labels=['Wi', 'Qwall', 'Hout-Hin', 'Qdiff'], autopct='%.2f')
# plt.title('Energy Balance Chart, Qin = ' + str(Qin) + ' kW')


plt.show()