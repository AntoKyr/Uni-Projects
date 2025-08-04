import numpy as np
from matplotlib import pyplot as plt
from tqdm import tqdm
from Engine_Model import EngineSim


Mode = 'test'
# can be map, point_Pint, point_phi_int or test

if Mode == 'map':

    Compute_Map = False
    # Points of op
    Nrange = np.arange(800, 6001, 200)
    Prange = np.arange(0.25, 1.01, 0.025)
    bmep_field = np.zeros((np.shape(Nrange)[0], np.shape(Prange)[0] + 2))
    bsfc_field = np.zeros((np.shape(Nrange)[0], np.shape(Prange)[0] + 2))

    # Calculating data-------------------------------------------------------------------------------------------------------
    if Compute_Map:
        Ncounter = 0
        for N in tqdm(Nrange, desc="Drawing map...", ascii=False, ncols=100):

            Pcounter = 0

            for Pint in tqdm(Prange, desc="Constant RPM sim...", ascii=False, ncols=100, leave=False):

                PerformanceParas = EngineSim(Pint=Pint, N=N)

                bmep_field[Ncounter, Pcounter] = PerformanceParas[0]
                bsfc_field[Ncounter, Pcounter] = PerformanceParas[-1]

                Pcounter += 1
            
            Ncounter += 1

        np.savez('MapData.npz', bmep_field=bmep_field, bsfc_field=bsfc_field)   # Save to avoid the need for multiple runs


    with np.load('MapData.npz', allow_pickle=True) as MapData:
        bmep_field = MapData['bmep_field']
        bsfc_field = MapData['bsfc_field']


    # Mapping data to grid points---------------------------------------------------------------------------------------------
    # Finding all evaluation points
    power_line = bmep_field[:,-3]
    increment = power_line + 10**-3
    maximum = np.ceil(max(power_line)+1)
    eval_points = np.sort(np.concatenate((np.arange(0, maximum+0.01, 0.1), power_line, increment)))
    bsfc_interp = np.zeros((np.shape(Nrange)[0], np.shape(eval_points)[0]))

    # Preparing data
    bmep_field[:,-2] = increment
    bmep_field[:,-1] = maximum

    # Interpolating Efficiencies
    for i in np.arange(0, np.shape(Nrange)[0]):
        # Eliminating non physical values
        bsfc_field[i,:][np.argwhere(bsfc_field[i,:]>10**6)] = 10**6
        bsfc_field[i,:][np.argwhere(bsfc_field[i,:]<0)] = 10**6

        # Interpolating
        bsfc_interp[i] = np.interp(x=eval_points, xp=bmep_field[i,:], fp=bsfc_field[i,:])

    plt.figure()
    plt.plot(bmep_field[0], bsfc_field[0], 'bo--', label='Interpolated Data')
    plt.plot(eval_points, bsfc_interp[0], 'r.--', label='Interpolation Points')
    plt.grid()
    plt.xlabel('bmep [bar]')
    plt.ylabel('bsfc [g/kWh]')
    plt.title('bsfc Interpolation at N = 800 rpm')
    plt.legend()
        
    effb_interp = np.transpose(bsfc_interp)

    # Plotting data----------------------------------------------------------------------------------------------------------
    contour_lvls = [165, 170, 175, 180, 190, 200, 210, 230, 250, 270, 300, 350]
    plt.figure()
    CS = plt.contour(Nrange, eval_points, effb_interp, cmap='jet', levels=contour_lvls, vmin=contour_lvls[0], vmax=contour_lvls[-1])
    plt.clabel(CS=CS, inline=True)
    plt.plot(Nrange, power_line, '-k', linewidth=5)
    plt.grid()
    plt.title('bsfc Contours [g/kWh]')
    plt.xlabel('Speed [rpm]')
    plt.ylabel('bmep [bar]')
    plt.show()



elif Mode == 'point_Pint':

    N = 2600
    bmep_target = 2.9

    bmep = 0
    Pint_1 = 0.97
    Pint_2 = 1
    accuracy = 10**-4
    
    while abs(bmep_target - bmep) > accuracy:
        
        bmep = EngineSim(Pint=(Pint_1 + Pint_2)/2, N=N, Updates=[True, False], int_phi_dur=101)[0]

        if bmep > bmep_target:
            Pint_2 = (Pint_1 + Pint_2)/2
        else:
            Pint_1 = (Pint_1 + Pint_2)/2

        Pint = (Pint_1 + Pint_2)/2
        print('Searching point progress: ' + str(100*accuracy/abs(bmep_target - bmep)) + '%')
    
    # Estimating Performance
    print('Pint = ' + str(Pint) + ' bar')
    PerformanceParas = EngineSim(Pint=Pint, N=N, Updates=[True, True], int_phi_dur=101)
    print('bmep = ' + str(PerformanceParas[0]) + ' bar')
    print('imep_net = ' + str(PerformanceParas[1]) + ' bar')
    print('imep_gross = ' + str(PerformanceParas[2]) + ' bar')
    print('pmep = ' + str(PerformanceParas[3]) + ' bar')
    print('effi = ' + str(PerformanceParas[4]))
    print('effb = ' + str(PerformanceParas[5]))
    print('ev = ' + str(PerformanceParas[6]))
    print('exh_f = ' + str(PerformanceParas[7]))
    print('dWidt = ' + str(PerformanceParas[8]) + ' kW')
    print('dWbdt = ' + str(PerformanceParas[9]) + ' kW')
    print('Ti = ' + str(PerformanceParas[10]) + ' Nm')
    print('Tb = ' + str(PerformanceParas[11]) + ' Nm')
    print('bsfc = ' + str(PerformanceParas[12]) + ' g/kWh')


# elif Mode == 'point_phi_int':
#     Pint = 1
#     N = 2600
#     bmep_target = 2.9

#     bmep = 0
#     int_phi_dur1 = 175
#     int_phi_dur2 = 180
#     accuracy = 10**-4
#     int_phi_dur=(int_phi_dur1 + int_phi_dur2)/2

#     while abs(bmep_target - bmep) > accuracy:
        
#         bmep = EngineSim(Pint=Pint, N=N, Updates=[True, False], int_phi_dur=int_phi_dur)[0]

#         if bmep > bmep_target:
#             int_phi_dur2 = (int_phi_dur1 + int_phi_dur2)/2
#         else:
#             int_phi_dur1 = (int_phi_dur1 + int_phi_dur2)/2

#         int_phi_dur = (int_phi_dur1 + int_phi_dur2)/2
#         print('Searching point progress: ' + str(100*accuracy/abs(bmep_target - bmep)) + '%')
    
#     # Estimating Performance
#     print('Pint = ' + str(Pint) + ' bar')
#     PerformanceParas = EngineSim(Pint=Pint, N=N, Updates=[True, True], int_phi_dur=int_phi_dur)
#     print('bmep = ' + str(PerformanceParas[0]) + ' bar')
#     print('imep_net = ' + str(PerformanceParas[1]) + ' bar')
#     print('imep_gross = ' + str(PerformanceParas[2]) + ' bar')
#     print('pmep = ' + str(PerformanceParas[3]) + ' bar')
#     print('effi = ' + str(PerformanceParas[4]))
#     print('effb = ' + str(PerformanceParas[5]))
#     print('ev = ' + str(PerformanceParas[6]))
#     print('exh_f = ' + str(PerformanceParas[7]))
#     print('dWidt = ' + str(PerformanceParas[8]) + ' kW')
#     print('dWbdt = ' + str(PerformanceParas[9]) + ' kW')
#     print('Ti = ' + str(PerformanceParas[10]) + ' Nm')
#     print('Tb = ' + str(PerformanceParas[11]) + ' Nm')
#     print('bsfc = ' + str(PerformanceParas[12]) + ' g/kWh')


        
elif Mode == 'test':

    Pint = 0.9782
    N = 2600
    PerformanceParas = EngineSim(Pint=Pint, N=N, Updates=[True, True], int_phi_dur=74)
    print('bmep = ' + str(PerformanceParas[0]) + ' bar')
    print('imep_net = ' + str(PerformanceParas[1]) + ' bar')
    print('imep_gross = ' + str(PerformanceParas[2]) + ' bar')
    print('pmep = ' + str(PerformanceParas[3]) + ' bar')
    print('effi = ' + str(PerformanceParas[4]))
    print('effb = ' + str(PerformanceParas[5]))
    print('ev = ' + str(PerformanceParas[6]))
    print('exh_f = ' + str(PerformanceParas[7]))
    print('dWidt = ' + str(PerformanceParas[8]) + ' kW')
    print('dWbdt = ' + str(PerformanceParas[9]) + ' kW')
    print('Ti = ' + str(PerformanceParas[10]) + ' Nm')
    print('Tb = ' + str(PerformanceParas[11]) + ' Nm')
    print('bsfc = ' + str(PerformanceParas[12]) + ' g/kWh')
