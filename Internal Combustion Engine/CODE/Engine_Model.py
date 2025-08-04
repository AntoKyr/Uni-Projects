import numpy as np
from matplotlib import pyplot as plt
from tqdm import tqdm


def EngineSim(Pint, N, int_phi_dur = 272, Updates=[False, False]):
    # Pint in bars, N in rpm
    ## Problem Parameters--------------------------------------------------------------------------------------------------------------------------------------------------------

    # UI Configuration
    Cycle_bar = not Updates[0]
    leave_bar = Updates[1]

    if Cycle_bar:
        Phase_bars = True
        leave_bar = False
    elif (not Cycle_bar) and (leave_bar):
        Phase_bars = False
    else:
        Phase_bars = True


    # Simulation Constants
    NoR = 10 # number of cycles

    # Operation Point
    Pint = 10**5 * Pint
    omega = N*360/60 #degrees/s
    lamda_nom = 1

    # Cylinder Geometry Constants
    b = 0.085 #m
    e = 11 
    Ap = (b/2)**2 * np.pi #m**2
    Ac = 1.15*Ap #m**2
    Acyl = Ap
    r = 0.045 #m
    l = 0.145 #m
    Vd = 2*r*Acyl #m**3
    Vh = e*Vd/(e-1) #m**3
    Vc = Vh - Vd #m**3
    ls = r/l 
    hskirt = 0.01 #m

    # Heat Transfer Constants
    cm = 2*r*omega/360 #m/s
    Tw = 400 #K

    # Combustion Constants
    cu = 2*cm
    a = 6.908
    Hu = 43.4 * 10**6 #J/kg
    comb_start = 693 #degress
    comb_acon = [0.596, 0.429, 1.355, 1.115, 0.964, 1.076, 1.046, 1.007]
    comb_bcon = [2.48, 0.031, -18.49, -0.346, 75.56, -2.534*10**-4, -4.075*10**-7, 0.004]
    comb_duration_A = 58
    comb_m_A = 2.3
    imep_A = 3.65114737292416 * 10**5
    exh_f_A = 0.20699323691014618
    comb_start_A = 715
    N_A = 2000

    # Justi Constants
    Tbez = 273.15 #K
    AF = 14.7

    # Gas Constant
    R = 287 # J/kgK

    # Valve Constants
    maxlift = 0.01 #m
    no_val = 4 / 2 

    int_phi_o = 365 #+ 80 #degrees
    int_phi_c = int_phi_o + int_phi_dur #+ 30 #degrees
    int_d = 0.035 #m
    int_vl_rd = np.linspace(0, 0.35, 15)
    int_Cf = np.array([0, 0.055, 0.124, 0.19, 0.253, 0.313, 0.365, 0.413, 0.453, 0.48, 0.503, 0.515, 0.521, 0.524, 0.525])

    exh_phi_o = 140 #+ 30 #degrees
    exh_phi_c = 360 #+ 110 #degrees
    exh_d = 0.03 #m
    exh_vl_rd = np.linspace(0, 0.49, 13)
    exh_Cf = np.array([0, 0.122, 0.279, 0.413, 0.518, 0.592, 0.63, 0.644, 0.652, 0.657, 0.66, 0.664, 0.664])

    # Estimating fmep
    fmep = np.polyval(p=np.polyfit(x=[1500, 2000, 2500, 3000, 4000, 4500, 5000], y=3*np.array([1.7, 1.9, 2.2, 2.55, 3.45, 3.9, 4.4])/4.65, deg=2), x=N)*10**5

    ## Models---------------------------------------------------------------------------------------------------------------------------------------------------------

    # Cylinder Volume Modeling
    def Cyl_Vol(theta):
        theta = np.pi * theta / 180
        Vol = Vc + Acyl*(r + l - np.sqrt(l**2 - r**2 * np.sin(theta)**2) - r*np.cos(theta))
        dV = omega * Acyl * r * (np.sin(theta) + 0.5*ls*np.sin(2*theta) / np.sqrt(1 - ls**2 * np.sin(2*theta)**2))
        return [Vol, dV]


    # Heat Transfer Modeling
    def Heat_Transfer(sys_status, theta, T, P, TPV=0, P0=0):

        if sys_status == 'closed':
            C1 = 2.28 + 0.308 * cu/cm
        elif sys_status == 'open':
            C1 = 6.18 + 0.417 * cu/cm
        
        TPV=0
        C2 = 0.00324
        theta = np.pi * theta / 180
        As = b * np.pi * (r + l - np.sqrt(l**2 - r**2 * np.sin(theta)**2) - r*np.cos(theta)) + b * np.pi * hskirt
        v = C1*cm + C2*Vh*TPV*(P - P0)
        a = 127.93 * b**(-0.2) * abs(P * v / 10**5)**0.8 * abs(T)**-0.53
        return [a*(Ap + Ac + As)*(Tw - T), a]


    # Combustion Modeling
    def Combustion(mf, comb_start, theta, comb_duration, omega):
        dtheta = theta - comb_start 
        dQbdt = omega * mf * Hu * a * (comb_m + 1) * (dtheta / comb_duration)**comb_m * np.exp(-a * (dtheta / comb_duration)**(comb_m + 1)) / comb_duration
        dmfdtheta = mf * a * (comb_m + 1) * (dtheta / comb_duration)**comb_m * np.exp(-a * (dtheta / comb_duration)**(comb_m + 1)) / comb_duration
        return [dQbdt, dmfdtheta]


    # Justi's Function
    def Justi(T, lamda):
        F = 1 / lamda
        u = 10**3 * 0.1445 * (1356.8 + (489.6 + 46.4 / lamda**0.93) * (T - Tbez) * 10**-2 + (7.768 + 3.36 / lamda**0.8) * (T - Tbez)**2 * 10**-4 - (0.0975 + 0.0485 / lamda**0.75) * (T - Tbez)**3 * 10**-6)
        dudT = 10**3 * 0.1445 * ((489.6 + 46.4 / lamda**0.93) * 10**-2 + 2 * (7.768 + 3.36 / lamda**0.8) * (T - Tbez) * 10**-4 - 3 * (0.0975 + 0.0485 / lamda**0.75) * (T - Tbez)**2 * 10**-6)
        dudF = 10**3 * 0.1445 * (46.4 * (T - Tbez) * 10**-2 * 0.93 * F**-0.07 + 3.36 * (T - Tbez)**2 * 10**-4 * 0.8 * F**-0.2 - 0.0485 * (T - Tbez)**3 *10**-6 * 0.75 * F**-0.25)
        return [u, dudT, dudF]


    # Valve Flow Modeling
    def Valve_Flow(T, P1, P2, Avalve, Cf, gamma):
        pp = P2/P1
        if pp < (2/(gamma + 1))**(gamma/(gamma-1)):
            return no_val * P1 / (R*T) * Cf * Avalve * np.sqrt(gamma * R * T) * (2/(gamma+1)) ** ((gamma+1)/(2*gamma-2)) # sonic flow
        else:
            return no_val * Cf * Avalve * np.sqrt(P1**2 / (R*T)) * np.sqrt( 2*gamma * ( pp ** (2/gamma) - pp ** ((gamma+1)/gamma) ) / (gamma-1) ) #subsonic flow


    # Valve Opening Function
    def Valve_Opening(valve_type, phi):

        if valve_type == 'int':
            valve_d = int_d
            phi_o = int_phi_o
            phi_c = int_phi_c
            vl_rd = int_vl_rd
            Cfs = int_Cf

        elif valve_type == 'exh':
            valve_d = exh_d
            phi_o = exh_phi_o
            phi_c = exh_phi_c
            vl_rd = exh_vl_rd
            Cfs = exh_Cf

        lift = 0.5 * maxlift * (1 - np.cos(2 * np.pi * (phi - phi_o) / (phi_c - phi_o)))
        Cf = np.interp(lift/maxlift,vl_rd,Cfs)
        
        if lift < valve_d/4:
            Avalve = np.pi * valve_d * lift
        else:
            Avalve = valve_d**2 * np.pi / 4
        
        return [Avalve, Cf]


    # Enthalpy Calculation
    def Enthalpy(T,lamda):
        u = 10**3 * 0.1445 * (1356.8 + (489.6 + 46.4 / lamda**0.93) * (T - Tbez) * 10**-2 + (7.768 + 3.36 / lamda**0.8) * (T - Tbez)**2 * 10**-4 - (0.0975 + 0.0485 / lamda**0.75) * (T - Tbez)**3 * 10**-6)
        return u + R*T


    # Gamma Calculation
    def Gamma(T, lamda):
        cv = 10**3 * 0.1445 * ((489.6 + 46.4 / lamda**0.93) * 10**-2 + 2 * (7.768 + 3.36 / lamda**0.8) * (T - Tbez) * 10**-4 - 3 * (0.0975 + 0.0485 / lamda**0.75) * (T - Tbez)**2 * 10**-6)
        cp = cv + R
        return cp/cv


    # Temperature Gradient
    def TempGrad(P, dVdt, dQwdt, dQbdt, dmdt, h, u, dmfdt, m, dudf, dFdt, cv, Overlap=False):
        if Overlap:
            return (-P * dVdt + (dQwdt + dQbdt + dmdt@h - u*(sum(dmdt) + dmfdt)) / m - dudf * dFdt) / cv
        else:
            return (-P * dVdt + (dQwdt + dQbdt + dmdt*h - u*(dmdt + dmfdt)) / m - dudf * dFdt) / cv


    # Combustion Constants Calculation
    def Comb_Constants(comb_start, exh_f, N, imep):

        if imep<0:
            imep = imep_A
        gzzp = (comb_acon[0]+comb_bcon[0]*comb_start**-0.5)/(comb_acon[0]+comb_bcon[0]*comb_start_A**-0.5)
        gxrg = (comb_acon[1]+comb_bcon[1]*exh_f)/(comb_acon[1]+comb_bcon[1]*exh_f_A)
        gn = (comb_acon[2]+comb_bcon[2]*N**-0.5)/(comb_acon[2]+comb_bcon[2]*N_A**-0.5)
        gwi = (comb_acon[3]+comb_bcon[3]*imep**-0.5)/(comb_acon[3]+comb_bcon[3]*imep_A**-0.5)
        hzzp = (comb_acon[4]+comb_bcon[4]*comb_start**-2)/(comb_acon[4]+comb_bcon[4]*comb_start_A**-2)
        hxrg = (comb_acon[5]+comb_bcon[5]*exh_f**2)/(comb_acon[5]+comb_bcon[5]*exh_f_A**2)
        hn = (comb_acon[6]+comb_bcon[6]*N**1.5)/(comb_acon[6]+comb_bcon[6]*N_A*1.5)
        hwi = (comb_acon[7]+comb_bcon[7]*np.log(imep))/(comb_acon[7]+comb_bcon[7]*np.log(imep_A))
        comb_duration = comb_duration_A * gzzp * gxrg * gn * gwi
        comb_m = comb_m_A * hzzp * hxrg * hn * hwi
        return [comb_duration, comb_m]



    ## Internal Combustion Engine Cylinder Simulation-----------------------------------------------------------------------------------------------------

    # Discretization
    disc_combustion = 5
    disc_closed = 2
    disc_open = 40

    # Function for creating discrete field
    def discrete_field(start, finnish, disc):
        return np.linspace(start, finnish-1/disc, (finnish - start)*disc)

    # Data Collection Bins
    imep_net = np.zeros(NoR)
    imep_gross = np.zeros(NoR)
    pmep = np.zeros(NoR)
    bmep = np.zeros(NoR)
    ev = np.zeros(NoR)
    exh_f = np.zeros(NoR)
    dWidt = np.zeros(NoR)
    dWbdt = np.zeros(NoR)
    dQindt = np.zeros(NoR)
    Ti = np.zeros(NoR)
    Tb = np.zeros(NoR)
    effi = np.zeros(NoR)
    effb = np.zeros(NoR)
    bsfc = np.zeros(NoR)
    comb_dur_data = np.zeros(NoR)
    comb_m_data = np.zeros(NoR)
    Energy_Fracts = np.zeros((NoR,4))

    # Start Conditions
    T = 1100
    TPV = 0
    T0 = T
    P = 1 * 10**5
    maircyl = 0
    mexhcyl = P*Vc / (R*T)
    mfcyl = 0
    mfinj = 0
    mexhint = 0
    mcyl = mexhcyl + maircyl
    lamda = (maircyl + mexhcyl) / (AF*mfcyl + mexhcyl)

    # Combustion Start Conditions
    comb_const = Comb_Constants(comb_start=comb_start, exh_f=exh_f_A, N=N, imep=imep_A)
    comb_duration = round(comb_const[0])
    comb_m = comb_const[1]

    # Boundary Conditions
    Tint = 300
    Pexh = 10**5
    IntAirEnthalpy = Enthalpy(T=Tint, lamda=10**8)



    for CR in tqdm(np.arange(0, NoR), desc="Simulating Engine...", ascii=False, ncols=100, disable=Cycle_bar, leave=leave_bar):

        # Setting up next cycle
        comb_end = comb_start + comb_duration - 720

        # Injection Constants
        inj_start = 540
        inj_end = inj_start + 5 #degrees

        compression_phis_1 = discrete_field(int_phi_c, inj_start, disc_closed)
        injection_phis = discrete_field(inj_start, inj_end, disc_closed)
        exhaust_phis = discrete_field(exh_phi_o, min(exh_phi_c, int_phi_o), disc_open)
        intake_phis = discrete_field(max(exh_phi_c, int_phi_o), int_phi_c, disc_open)
        combustion_phis_2 = discrete_field(comb_start, 720, disc_combustion)
        combustion_phis_1 = discrete_field(0, comb_end, disc_combustion)
        power_phis_1 = discrete_field(comb_end, exh_phi_o, disc_closed)

        if inj_end < comb_start:
            compression_phis_2 = discrete_field(inj_end, comb_start, disc_closed)
        else:
            compression_phis_2 = []

        if exh_phi_c < int_phi_o:
            power_phis_2 = discrete_field(exh_phi_c, int_phi_o, disc_closed)
            overlap_phis = []
        elif exh_phi_c == int_phi_o:
            power_phis_2 = []
            overlap_phis = []
        elif exh_phi_c > int_phi_o:
            power_phis_2 = []
            overlap_phis = discrete_field(int_phi_o, exh_phi_c, disc_open)

        phis = np.concatenate((combustion_phis_1, power_phis_1, exhaust_phis, power_phis_2, overlap_phis, intake_phis, compression_phis_1, injection_phis, compression_phis_2, combustion_phis_2))
        field_size = np.shape(phis)[0]

        # Surprise Tool that will help us later
        pmep_indx= np.array([np.argwhere(phis>180)[0,0], np.argwhere(phis>540)[0,0]])

        # More data collection
        P_data = np.zeros(field_size)
        T_data = np.zeros(field_size)
        V_data = np.zeros(field_size)
        HT_data = np.zeros(field_size)
        HTC_data = np.zeros(field_size)
        Comb_data = np.zeros(field_size)
        mair_data = np.zeros(field_size)
        mexh_data = np.zeros(field_size)
        mf_data = np.zeros(field_size)
        int_mflow_data = np.zeros(field_size)
        exh_mflow_data = np.zeros(field_size)
        exhint_data = np.zeros(field_size)
        AV_data = np.zeros(field_size)
        Cf_data = np.zeros(field_size)

        counter = 0
        Hdiff = 0
        Qw = 0



        # Combustion
        dphi = 1/disc_combustion
        for phi in tqdm(combustion_phis_1, desc="Combustion", ascii=False, ncols=100, disable=Phase_bars, leave=False): 
            
            # Computing parallel step
            Volumes = Cyl_Vol(phi)
            P0 = mcyl * R * T0 / Volumes[0]
            Justi_data0 = Justi(T=T0, lamda=lamda)
            dQwdt0 = Heat_Transfer(sys_status='closed', theta=phi, T=T0, P=P0)[0]
            dTdt0 = TempGrad(P=P0, dVdt=Volumes[1], dQwdt=dQwdt0, dQbdt=0, dmdt=0, h=0, u=Justi_data0[0], dmfdt=0, m=mcyl, dudf=Justi_data0[2], dFdt=0, cv=Justi_data0[1])
            T0 = T0 + dphi*dTdt0/omega

            # Computing current step
            P = mcyl * R * T / Volumes[0]
            Justi_data = Justi(T=T, lamda=lamda)
            Woschni = Heat_Transfer(sys_status='closed', theta=phi, T=T, P=P, TPV=TPV, P0=P0)
            dQwdt = Woschni[0]
            Comb_rates = Combustion(mf=mfinj, comb_start=comb_start-720, theta=phi, comb_duration=comb_duration, omega=omega)
            dTdt = TempGrad(P=P, dVdt=Volumes[1], dQwdt=dQwdt, dQbdt=Comb_rates[0], dmdt=0, h=0, u=Justi_data[0], dmfdt=0, m=mcyl, dudf=Justi_data[2], dFdt=0, cv=Justi_data[1])

            # Saving data
            P_data[counter] = P
            T_data[counter] = T
            V_data[counter] = Volumes[0]
            HT_data[counter] = dQwdt
            HTC_data[counter] = Woschni[1]
            Comb_data[counter] = Comb_rates[0]
            mair_data[counter] = maircyl
            mexh_data[counter] = mexhcyl
            mf_data[counter] = mfcyl
            int_mflow_data[counter] = 0
            exh_mflow_data[counter] = 0
            exhint_data[counter] = mexhint
            AV_data[counter] = 0
            Cf_data[counter] = 0

            # Computing next step
            maircyl = maircyl - AF * Comb_rates[1] * dphi
            mexhcyl = mexhcyl + AF * Comb_rates[1] * dphi
            mfcyl = mfcyl - Comb_rates[1] * dphi
            T = T + dTdt * dphi/omega
            Qw = Qw - dQwdt * dphi/omega
            counter += 1
        


        # Approximating that the whole of the fuel mass has burnt
        mexhcyl = mexhcyl + (AF + 1) * mfcyl
        maircyl = maircyl - AF * mfcyl
        mfcyl = 0
        yair = maircyl / mcyl
        yexh = mexhcyl / mcyl
        


        # Power
        dphi = 1/disc_closed
        for phi in tqdm(power_phis_1, desc="Power_1", ascii=False, ncols=100, disable=Phase_bars, leave=False): 

            # Computing current step
            Volumes = Cyl_Vol(phi)
            P = mcyl * R * T / Volumes[0]
            Justi_data = Justi(T=T, lamda=lamda)
            Woschni = Heat_Transfer(sys_status='closed', theta=phi, T=T, P=P)
            dQwdt = Woschni[0]
            dTdt = TempGrad(P=P, dVdt=Volumes[1], dQwdt=dQwdt, dQbdt=0, dmdt=0, h=0, u=Justi_data[0], dmfdt=0, m=mcyl, dudf=Justi_data[2], dFdt=0, cv=Justi_data[1])

            # Saving data
            P_data[counter] = P
            T_data[counter] = T
            V_data[counter] = Volumes[0]
            HT_data[counter] = dQwdt
            HTC_data[counter] = Woschni[1]
            Comb_data[counter] = 0
            mair_data[counter] = maircyl
            mexh_data[counter] = mexhcyl
            mf_data[counter] = mfcyl
            int_mflow_data[counter] = 0
            exh_mflow_data[counter] = 0
            exhint_data[counter] = mexhint
            AV_data[counter] = 0
            Cf_data[counter] = 0

            # Computing next step    
            T = T + dTdt * dphi/omega
            Qw = Qw - dQwdt * dphi/omega
            counter += 1    



        # Exhaust 
        dphi = 1/disc_open
        for phi in tqdm(exhaust_phis, desc="Exhaust_2", ascii=False, ncols=100, disable=Phase_bars, leave=False): 

            # Computing current step
            mcyl = maircyl + mexhcyl
            Volumes = Cyl_Vol(phi)
            P = mcyl * R * T / Volumes[0]
            Justi_data = Justi(T=T, lamda=lamda)
            Valve_data = Valve_Opening(valve_type='exh', phi=phi)
            Woschni = Heat_Transfer(sys_status='open', theta=phi, T=T, P=P)
            dQwdt = Woschni[0]
            
            if P/Pexh >= 1:
                dmdt = -Valve_Flow(T=T, P1=P, P2=Pexh, Avalve=Valve_data[0], Cf=Valve_data[1], gamma=Gamma(T=T, lamda=lamda))
                h = Enthalpy(T=T, lamda=lamda)

            else:
                dmdt = Valve_Flow(T=T, P1=Pexh, P2=P, Avalve=Valve_data[0], Cf=Valve_data[1], gamma=Gamma(T=T, lamda=lamda))
                h = Enthalpy(T=T, lamda=lamda)

            dTdt = TempGrad(P=P, dVdt=Volumes[1], dQwdt=dQwdt, dQbdt=0, dmdt=dmdt, h=h, u=Justi_data[0], dmfdt=0, m=mcyl, dudf=Justi_data[2], dFdt=0, cv=Justi_data[1])

            # Saving data
            P_data[counter] = P
            T_data[counter] = T
            V_data[counter] = Volumes[0]
            HT_data[counter] = dQwdt
            HTC_data[counter] = Woschni[1]
            Comb_data[counter] = 0
            mair_data[counter] = maircyl
            mexh_data[counter] = mexhcyl
            mf_data[counter] = mfcyl
            int_mflow_data[counter] = 0
            exh_mflow_data[counter] = dmdt
            exhint_data[counter] = mexhint
            AV_data[counter] = Valve_data[0]
            Cf_data[counter] = Valve_data[1]

            # Computing next step
            maircyl = maircyl + yair * dmdt * dphi/omega
            mexhcyl = mexhcyl + yexh * dmdt * dphi/omega
            T = T + dTdt * dphi/omega
            Qw = Qw - dQwdt * dphi/omega
            Hdiff = Hdiff + h*dmdt * dphi/omega
            counter += 1

        

        # Valve Overlap
        dphi = 1/disc_open
        for phi in tqdm(overlap_phis, desc="Overlap", ascii=False, ncols=100, disable=Phase_bars, leave=False): 

            # Computing current step       
            mcyl = maircyl + mexhcyl
            Volumes = Cyl_Vol(phi)
            P = mcyl * R * T / Volumes[0]
            Justi_data = Justi(T=T, lamda=lamda)
            Woschni = Heat_Transfer(sys_status='open', theta=phi, T=T, P=P)
            dQwdt = Woschni[0]
            Valve_data = Valve_Opening(valve_type='exh', phi=phi)

            if P/Pexh >= 1:
                dmdte = -Valve_Flow(T=T, P1=P, P2=Pexh, Avalve=Valve_data[0], Cf=Valve_data[1], gamma=Gamma(T=T, lamda=lamda))
                he = Enthalpy(T=T, lamda=lamda)

            else:
                dmdte = Valve_Flow(T=T, P1=Pexh, P2=P, Avalve=Valve_data[0], Cf=Valve_data[1], gamma=Gamma(T=T, lamda=lamda_nom))
                he = Enthalpy(T=T, lamda=lamda_nom)
            
            Valve_data = Valve_Opening(valve_type='int', phi=phi)

            if P/Pint < 1:

                if mexhint <= 0:
                    mexhint = 0
                    yexh = 0
                    dmdti = Valve_Flow(T=Tint, P1=Pint, P2=P, Avalve=Valve_data[0], Cf=Valve_data[1], gamma=1.4)
                    hi = IntAirEnthalpy
                    dFdt = (dmdte * (mexhcyl + maircyl) - mexhcyl * (dmdte + dmdti)) / (mexhcyl + maircyl)**2 

                else:
                    yexh = 1/lamda_nom
                    dmdti = Valve_Flow(T=T, P1=Pint, P2=P, Avalve=Valve_data[0], Cf=Valve_data[1], gamma=Gamma(T=T, lamda=lamda_nom))
                    hi = Enthalpy(T=T, lamda=lamda_nom)
                    dFdt = 0

            else:
                yexh = mexhcyl/(maircyl + mexhcyl)
                dmdti = - Valve_Flow(T=T, P1=P, P2=Pint, Avalve=Valve_data[0], Cf=Valve_data[1], gamma=Gamma(T=T, lamda=1/yexh))
                hi = Enthalpy(T=T, lamda=1/yexh)
                dFdt = 0
            

            dmdt = np.array([dmdti, dmdte])
            h = np.array([hi, he])
            dTdt = TempGrad(P=P, dVdt=Volumes[1], dQwdt=dQwdt, dQbdt=0, dmdt=dmdt, h=h, u=Justi_data[0], dmfdt=0, m=mcyl, dudf=Justi_data[2], dFdt=dFdt, cv=Justi_data[1], Overlap=True)

            # Saving data
            P_data[counter] = P
            T_data[counter] = T
            V_data[counter] = Volumes[0]
            HT_data[counter] = dQwdt
            HTC_data[counter] = Woschni[1]
            Comb_data[counter] = 0
            mair_data[counter] = maircyl
            mexh_data[counter] = mexhcyl
            mf_data[counter] = mfcyl
            int_mflow_data[counter] = dmdti
            exh_mflow_data[counter] = dmdte
            exhint_data[counter] = mexhint
            AV_data[counter] = Valve_data[0]
            Cf_data[counter] = Valve_data[1]

            # Computing next step
            yair = 1 - yexh        # Its important to note that at this stage yexh and yair refer to the gas in the intake
            mexhint = mexhint - yexh*dmdti * dphi/omega
            maircyl = maircyl + (yair*dmdti + (1-1/lamda)*dmdte) * dphi/omega
            mexhcyl = mexhcyl + (yexh*dmdti + (1/lamda)*dmdte) * dphi/omega
            lamda = (mexhcyl + maircyl) / mexhcyl
            T = T + dTdt * dphi/omega
            Qw = Qw - dQwdt * dphi/omega
            Hdiff = Hdiff + (he*dmdte + hi*dmdti) * dphi/omega
            counter += 1

            

        # Power
        dphi = 1/disc_closed
        for phi in tqdm(power_phis_2, desc="Power_1", ascii=False, ncols=100, disable=Phase_bars, leave=False): 

            # Computing current step
            Volumes = Cyl_Vol(phi)
            P = mcyl * R * T / Volumes[0]
            Justi_data = Justi(T=T, lamda=lamda)
            Woschni = Heat_Transfer(sys_status='closed', theta=phi, T=T, P=P)
            dQwdt = Woschni[0]
            dTdt = TempGrad(P=P, dVdt=Volumes[1], dQwdt=dQwdt, dQbdt=0, dmdt=0, h=0, u=Justi_data[0], dmfdt=0, m=mcyl, dudf=Justi_data[2], dFdt=0, cv=Justi_data[1])

            # Saving data
            P_data[counter] = P
            T_data[counter] = T
            V_data[counter] = Volumes[0]
            HT_data[counter] = dQwdt
            HTC_data[counter] = Woschni[1]
            Comb_data[counter] = 0
            mair_data[counter] = maircyl
            mexh_data[counter] = mexhcyl
            mf_data[counter] = mfcyl
            int_mflow_data[counter] = 0
            exh_mflow_data[counter] = 0
            exhint_data[counter] = mexhint
            AV_data[counter] = 0
            Cf_data[counter] = 0

            # Computing next step    
            T = T + dTdt * dphi/omega
            Qw = Qw - dQwdt * dphi/omega
            counter += 1


        
        # Intake
        dphi = 1/disc_open
        for phi in tqdm(intake_phis, desc="Intake", ascii=False, ncols=100, disable=Phase_bars, leave=False): 

            # Computing current step
            mcyl = maircyl + mexhcyl
            lamda = (mexhcyl + maircyl) / mexhcyl 
            Volumes = Cyl_Vol(phi)
            P = mcyl * R * T / Volumes[0]
            Justi_data = Justi(T=T, lamda=lamda)
            Valve_data = Valve_Opening(valve_type='int', phi=phi)
            Woschni = Heat_Transfer(sys_status='open', theta=phi, T=T, P=P)
            dQwdt = Woschni[0]
            
            if P/Pint < 1:
                
                if mexhint <= 0:
                    mexhint = 0
                    yexh = 0
                    dmdt = Valve_Flow(T=Tint, P1=Pint, P2=P, Avalve=Valve_data[0], Cf=Valve_data[1], gamma=1.4)
                    h = IntAirEnthalpy
                    dFdt = -mexhcyl * dmdt / (mexhcyl + maircyl)**2

                else:
                    yexh = 1/lamda_nom
                    dmdt = Valve_Flow(T=T, P1=Pint, P2=P, Avalve=Valve_data[0], Cf=Valve_data[1], gamma=Gamma(T=T, lamda=lamda_nom))
                    h = Enthalpy(T=T, lamda=lamda_nom)
                    dFdt = 0

            else:
                yexh = mexhcyl/(maircyl + mexhcyl)
                dmdt = - Valve_Flow(T=T, P1=P, P2=Pint, Avalve=Valve_data[0], Cf=Valve_data[1], gamma=Gamma(T=T, lamda=1/yexh))
                h = Enthalpy(T=T, lamda=1/yexh)
                dFdt = 0
            
            dTdt = TempGrad(P=P, dVdt=Volumes[1], dQwdt=dQwdt, dQbdt=0, dmdt=dmdt, h=h, u=Justi_data[0], dmfdt=0, m=mcyl, dudf=Justi_data[2], dFdt=dFdt, cv=Justi_data[1])

            # Saving data
            P_data[counter] = P
            T_data[counter] = T
            V_data[counter] = Volumes[0]
            HT_data[counter] = dQwdt
            HTC_data[counter] = Woschni[1]
            Comb_data[counter] = 0
            mair_data[counter] = maircyl
            mexh_data[counter] = mexhcyl
            mf_data[counter] = mfcyl
            int_mflow_data[counter] = dmdt
            exh_mflow_data[counter] = 0
            exhint_data[counter] = mexhint
            AV_data[counter] = Valve_data[0]
            Cf_data[counter] = Valve_data[1]

            # Computing next step
            yair = 1 - yexh
            mexhint = mexhint - yexh*dmdt * dphi/omega
            maircyl = maircyl + yair*dmdt * dphi/omega
            mexhcyl = mexhcyl + yexh*dmdt * dphi/omega
            T = T + dTdt * dphi/omega
            Qw = Qw - dQwdt * dphi/omega
            Hdiff = Hdiff + h*dmdt * dphi/omega
            counter += 1

        

        # Computing Heat Transfer Parameter for Combustion Phase
        TPV = T / (P*Volumes[0])
        # Exhaust Fraction and Filling Efficiency
        ev[CR] = maircyl * R * Tint / (Pint * Vd)
        exh_f[CR] = mexhcyl / (mexhcyl + maircyl)



        # Compression
        dphi = 1/disc_closed
        for phi in tqdm(compression_phis_1, desc="Compression_1", ascii=False, ncols=100, disable=Phase_bars, leave=False): 
            
            # Computing current step
            Volumes = Cyl_Vol(phi)
            P = mcyl * R * T / Volumes[0]
            Justi_data = Justi(T=T, lamda=lamda)
            Woschni = Heat_Transfer(sys_status='closed', theta=phi, T=T, P=P)
            dQwdt = Woschni[0]
            dTdt = TempGrad(P=P, dVdt=Volumes[1], dQwdt=dQwdt, dQbdt=0, dmdt=0, h=0, u=Justi_data[0], dmfdt=0, m=mcyl, dudf=Justi_data[2], dFdt=0, cv=Justi_data[1])

            # Saving data
            P_data[counter] = P
            T_data[counter] = T
            V_data[counter] = Volumes[0]
            HT_data[counter] = dQwdt
            HTC_data[counter] = Woschni[1]
            Comb_data[counter] = 0
            mair_data[counter] = maircyl
            mexh_data[counter] = mexhcyl
            mf_data[counter] = mfcyl
            int_mflow_data[counter] = 0
            exh_mflow_data[counter] = 0
            exhint_data[counter] = mexhint
            AV_data[counter] = 0
            Cf_data[counter] = 0

            # Computing next step
            T = T + dTdt * dphi/omega
            Qw = Qw - dQwdt * dphi/omega
            counter += 1



        # Computing Injection Parameters
        mfinj = maircyl / (AF*lamda_nom)
        inj_duration = inj_end - inj_start
        inj_rate = mfinj / (inj_duration/omega)
        Qin = mfinj * Hu



        # Injection
        dphi = 1/disc_closed
        for phi in tqdm(injection_phis, desc="Injection", ascii=False, ncols=100, disable=Phase_bars, leave=False): 
            
            # Computing current step
            mcyl = maircyl + mexhcyl + mfcyl
            Volumes = Cyl_Vol(phi)
            P = mcyl * R * T / Volumes[0]
            Justi_data = Justi(T=T, lamda=lamda)
            Woschni = Heat_Transfer(sys_status='closed', theta=phi, T=T, P=P)
            dQwdt = Woschni[0]
            dFdt = AF * inj_rate / maircyl
            dTdt = TempGrad(P=P, dVdt=Volumes[1], dQwdt=dQwdt, dQbdt=0, dmdt=0, h=0, u=Justi_data[0], dmfdt=inj_rate, m=mcyl, dudf=Justi_data[2], dFdt=dFdt, cv=Justi_data[1])

            # Saving data
            P_data[counter] = P
            T_data[counter] = T
            V_data[counter] = Volumes[0]
            HT_data[counter] = dQwdt
            HTC_data[counter] = Woschni[1]
            Comb_data[counter] = 0
            mair_data[counter] = maircyl
            mexh_data[counter] = mexhcyl
            mf_data[counter] = mfcyl
            int_mflow_data[counter] = 0
            exh_mflow_data[counter] = 0
            exhint_data[counter] = mexhint
            AV_data[counter] = 0
            Cf_data[counter] = 0

            # Computing next step
            mfcyl = mfcyl + inj_rate * dphi/omega
            T = T + dTdt * dphi/omega
            Qw = Qw - dQwdt * dphi/omega
            lamda = (maircyl + mexhcyl) / (AF*mfcyl + mexhcyl)
            counter += 1



        # Compression
        dphi = 1/disc_closed
        for phi in tqdm(compression_phis_2, desc="Compression_2", ascii=False, ncols=100, disable=Phase_bars, leave=False): 
            
            # Computing current step
            Volumes = Cyl_Vol(phi)
            P = mcyl * R * T / Volumes[0]
            Justi_data = Justi(T=T, lamda=lamda)
            Woschni = Heat_Transfer(sys_status='closed', theta=phi, T=T, P=P)
            dQwdt = Woschni[0]
            dTdt = TempGrad(P=P, dVdt=Volumes[1], dQwdt=dQwdt, dQbdt=0, dmdt=0, h=0, u=Justi_data[0], dmfdt=0, m=mcyl, dudf=Justi_data[2], dFdt=0, cv=Justi_data[1])

            # Saving data
            P_data[counter] = P
            T_data[counter] = T
            V_data[counter] = Volumes[0]
            HT_data[counter] = dQwdt
            HTC_data[counter] = Woschni[1]
            Comb_data[counter] = 0
            mair_data[counter] = maircyl
            mexh_data[counter] = mexhcyl
            mf_data[counter] = mfcyl
            int_mflow_data[counter] = 0
            exh_mflow_data[counter] = 0
            exhint_data[counter] = mexhint
            AV_data[counter] = 0
            Cf_data[counter] = 0

            # Computing next step
            T = T + dTdt * dphi/omega
            Qw = Qw - dQwdt * dphi/omega
            counter += 1



        # Initiating parallel simulation
        T0 = T



        # Combustion
        dphi = 1/disc_combustion
        for phi in tqdm(combustion_phis_2, desc="Combustion", ascii=False, ncols=100, disable=Phase_bars, leave=False): 
            
            # Computing parallel step
            Volumes = Cyl_Vol(phi)
            P0 = mcyl * R * T0 / Volumes[0]
            Justi_data0 = Justi(T=T0, lamda=lamda)
            dQwdt0 = Heat_Transfer(sys_status='closed', theta=phi, T=T0, P=P0)[0]
            dTdt0 = TempGrad(P=P0, dVdt=Volumes[1], dQwdt=dQwdt0, dQbdt=0, dmdt=0, h=0, u=Justi_data0[0], dmfdt=0, m=mcyl, dudf=Justi_data0[2], dFdt=0, cv=Justi_data0[1])
            T0 = T0 + dphi*dTdt0/omega

            # Computing current step
            P = mcyl * R * T / Volumes[0]
            Justi_data = Justi(T=T, lamda=lamda)
            Woschni = Heat_Transfer(sys_status='closed', theta=phi, T=T, P=P, TPV=TPV, P0=P0)
            dQwdt = Woschni[0]
            Comb_rates = Combustion(mf=mfinj, comb_start=comb_start, theta=phi, comb_duration=comb_duration, omega=omega)
            dTdt = TempGrad(P=P, dVdt=Volumes[1], dQwdt=dQwdt, dQbdt=Comb_rates[0], dmdt=0, h=0, u=Justi_data[0], dmfdt=0, m=mcyl, dudf=Justi_data[2], dFdt=0, cv=Justi_data[1])

            # Saving data
            P_data[counter] = P
            T_data[counter] = T
            V_data[counter] = Volumes[0]
            HT_data[counter] = dQwdt
            HTC_data[counter] = Woschni[1]
            Comb_data[counter] = Comb_rates[0]
            mair_data[counter] = maircyl
            mexh_data[counter] = mexhcyl
            mf_data[counter] = mfcyl
            int_mflow_data[counter] = 0
            exh_mflow_data[counter] = 0
            exhint_data[counter] = mexhint
            AV_data[counter] = 0
            Cf_data[counter] = 0

            # Computing next step
            maircyl = maircyl - AF * Comb_rates[1] * dphi
            mexhcyl = mexhcyl + AF * Comb_rates[1] * dphi
            mfcyl = mfcyl - Comb_rates[1] * dphi
            T = T + dTdt * dphi/omega
            Qw = Qw - dQwdt * dphi/omega
            counter += 1



        # Calculating performace indicators
        W1 = np.trapz(P_data[0:pmep_indx[0]+1], V_data[0:pmep_indx[0]+1])
        W2 = np.trapz(np.concatenate((P_data[pmep_indx[1]:], [P_data[0]])), np.concatenate((V_data[pmep_indx[1]:], [V_data[0]])))
        Wgross = W1 + W2
        Wnet = np.trapz(np.concatenate((P_data, [P_data[0]])), np.concatenate((V_data, [V_data[0]])))
        Wp = np.trapz(P_data[pmep_indx[0] : pmep_indx[1]+1], V_data[pmep_indx[0] : pmep_indx[1]+1])
        
        imep_net[CR] = Wnet / Vd
        imep_gross[CR] = Wgross / Vd
        pmep[CR] = Wp / Vd
        bmep[CR] = imep_gross[CR] - fmep
        dWidt[CR] = Wnet * omega/720
        dQindt[CR] = Qin * omega/720
        Energy_Fracts[CR,:] = np.array([Qin, Qw, Wnet, Hdiff]) * omega/720
        dWbdt[CR] = Vd * (N/60) * bmep[CR] / 2
        Tb[CR] = Vd * bmep[CR] / (4*np.pi)
        Ti[CR] = Vd * imep_net[CR] / (4*np.pi)
        effi[CR] = Wnet / Qin
        effb[CR] = dWbdt[CR] / dQindt[CR]
        bsfc[CR] = 3600 * 10**6 / (Hu * effb[CR])

        # Calculating new combustion constants
        comb_const = Comb_Constants(comb_start=comb_start, exh_f=exh_f[CR], N=N, imep=imep_net[CR])
        comb_duration = round(comb_const[0])
        comb_m = comb_const[1]
        comb_dur_data[CR] = comb_duration
        comb_m_data[CR] = comb_m


    ## Saving Data-------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Timings = [[int_phi_o, int_phi_c], [exh_phi_o, exh_phi_c], [comb_start, comb_end], [inj_start, inj_end]]
    np.savez('SimData.npz', NoR=NoR, phis=phis, Timings=Timings, P_data=P_data, T_data=T_data, HT_data=HT_data, HTC_data=HTC_data, Comb_data=Comb_data, mexh_data=mexh_data, 
            mair_data=mair_data, mf_data=mf_data, int_mflow_data=int_mflow_data, exh_mflow_data=exh_mflow_data, exhint_data=exhint_data, V_data=V_data,
            AV_data=AV_data, Cf_data=Cf_data, ev=ev, exh_f=exh_f, imep_gross=imep_gross, imep_net=imep_net, pmep=pmep, bmep=bmep, effi=effi, effb=effb,
            dWidt=dWidt, dWbdt=dWbdt, dQindt=dQindt, Ti=Ti, Tb=Tb, bsfc=bsfc, comb_dur_data=comb_dur_data, comb_m_data=comb_m_data, Energy_Fracts=Energy_Fracts)

    return [bmep[-1]/10**5 ,imep_net[-1]/10**5, imep_gross[-1]/10**5, -pmep[-1]/10**5, effi[-1], effb[-1], ev[-1], exh_f[-1], dWidt[-1]/10**3, dWbdt[-1]/10**3, Ti[-1], Tb[-1], bsfc[-1]]