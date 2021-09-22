'''
Conceptual Aircraft Design Tool
(for PRJ-22 and AP-701 courses)

Cap. Eng. Ney Rafael Secco (ney@ita.br)
Aircraft Design Department
Aeronautics Institute of Technology

07-2021

The code uses several historical regression from
aircraft design books to make a quick initial
sizing procedure.

Generally, the user should call only the 'analyze'
function from this module.
'''

# IMPORTS
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# CONSTANTS
ft2m = 0.3048
kt2ms = 0.514444
lb2N = 4.44822
nm2m = 1852.0
gravity = 9.81

#========================================
# MAIN FUNCTION

def analyze(airplane = None,
            print_log = False, # Plot results on the terminal screen
            plot = False, # Generate 3D plot of the aircraft
            W0_guess = None, # Guess for MTOW [N]
            T0_guess = None, # Guess for Takeoff total thrust [N]
            ):
    '''
    This is the main function that should be used for aircraft analysis.
    '''

    # Load standard airplane if none is provided
    if airplane is None:
        airplane = standard_airplane()

    # Use an average wing loading for transports
    # to estime W0_guess and T0_guess if none are provided
    if W0_guess is None:
        W0_guess = 5e3*airplane['S_w']
    if T0_guess is None:
        T0_guess = 0.3*W0_guess

    ### ADD CODE FROM SECTION 4.3 HERE ###


    if print_log:
        print('We', airplane['We']/gravity)
        print('Wf', airplane['Wf']/gravity)
        print('W0', airplane['W0']/gravity)
        print('T0', airplane['T0']/gravity)
        print('T0/W0', airplane['T0']/airplane['W0'])
        print('W0/S', airplane['W0']/airplane['S_w'])
        print('deltaS_wlan', airplane['deltaS_wlan'])
        print('xcg_fwd', airplane['xcg_fwd'])
        print('xcg_aft', airplane['xcg_aft'])
        print('xnp', airplane['xnp'])
        print('SM_fwd', airplane['SM_fwd'])
        print('SM_aft', airplane['SM_aft'])
        print('b_tank_b_w', airplane['b_tank_b_w'])
        print('CLv', airplane['CLv'])
        print('frac_nlg_fwd', airplane['frac_nlg_fwd'])
        print('frac_nlg_aft', airplane['frac_nlg_aft'])
        print('alpha_tipback [deg]', airplane['alpha_tipback']*180.0/np.pi)
        print('alpha_tailstrike [deg]', airplane['alpha_tailstrike']*180.0/np.pi)
        print('phi_overturn [deg]', airplane['phi_overturn']*180.0/np.pi)

    if plot:
        plot3d(airplane)

    return airplane

#========================================
# DISCIPLINE MODULES

def geometry(airplane):

    # Unpack dictionary
    S_w = airplane['S_w']
    AR_w = airplane['AR_w']
    taper_w = airplane['taper_w']
    sweep_w = airplane['sweep_w']
    dihedral_w = airplane['dihedral_w']
    xr_w = airplane['xr_w']
    zr_w = airplane['zr_w']
    Cht = airplane['Cht']
    AR_h = airplane['AR_h']
    taper_h = airplane['taper_h']
    sweep_h = airplane['sweep_h']
    dihedral_h = airplane['dihedral_h']
    Lc_h = airplane['Lc_h']
    zr_h = airplane['zr_h']
    Cvt = airplane['Cvt']
    AR_v = airplane['AR_v']
    taper_v = airplane['taper_v']
    sweep_v = airplane['sweep_v']
    Lb_v = airplane['Lb_v']
    zr_v = airplane['zr_v']

    ### Sizing of the wing:
    
    b_w = np.sqrt(AR_w*S_w)
    cr_w = 2*S_w/(b_w * (1 + taper_w))
    print(cr_w)
    ct_w = taper_w * cr_w
    yt_w = b_w/2
    xt_w = xr_w + yt_w * np.tan(sweep_w) + (cr_w - ct_w)/4
    zt_w = zr_w + yt_w * np.tan(dihedral_w)
    cm_w = (2*cr_w/3)*(1 + taper_w + taper_w**2)/(1 + taper_w)
    ym_w = (b_w/6)*(1 + 2*taper_w)/(1 + taper_w)
    xm_w = xr_w + ym_w * np.tan(sweep_w) + (cr_w - cm_w)/4
    zm_w = zr_w + ym_w * np.tan(dihedral_w)
    
    ### Sizing of the horizontal tail:
    
    L_h = Lc_h * cm_w
    S_h = S_w * Cht / Lc_h
    b_h = np.sqrt(AR_h * S_h)
    cr_h = 2 * S_h / (b_h * (1 + taper_h))
    ct_h = taper_h * cr_h
    cm_h = (2 * cr_h / 3) * (1 + taper_h + taper_h**2)/(1 + taper_h)
    xm_h = xm_w + L_h + (cm_w - cm_h)/4
    ym_h = (b_h/6) * (1 + 2 * taper_h) / (1 + taper_h)
    zm_h = zr_h + ym_h * np.tan(dihedral_h)
    xr_h = xm_h - ym_h * np.tan(sweep_h) + (cm_h - cr_h)/4
    yt_h = b_h/2
    xt_h = xr_h + yt_h * np.tan(sweep_h) + (cr_h - ct_h)/4
    zt_h = zr_h + yt_h * np.tan(dihedral_h)
    
    ### Sizing of the vertical tail:
    
    L_v = Lb_v * b_w
    S_v = S_w * Cvt / Lb_v
    b_v = np.sqrt(AR_v * S_v)
    cr_v = 2 * S_v / (b_v * (1 + taper_v))
    ct_v = taper_v * cr_v
    cm_v = (2*cr_v/3) * (1 + taper_v + taper_v**2) / (1 + taper_v)
    xm_v = xm_w + L_v + (cm_w - cm_v)/4
    zm_v = zr_v + b_v/3 + (1 + 2*taper_v)/(1 + taper_v)
    xr_v = xm_v - (zm_v - zr_v) * np.tan(sweep_v) + (cm_v - cr_v)/4
    zt_v = zr_v + b_v
    xt_v = xr_v + (zt_v - zr_v) * np.tan(sweep_v) + (cr_v - ct_v)/4

    # Update dictionary with new results
    airplane['b_w'] = b_w
    airplane['cr_w'] = cr_w
    airplane['xt_w'] = xt_w
    airplane['yt_w'] = yt_w
    airplane['zt_w'] = zt_w
    airplane['ct_w'] = ct_w
    airplane['xm_w'] = xm_w
    airplane['ym_w'] = ym_w
    airplane['zm_w'] = zm_w
    airplane['cm_w'] = cm_w
    airplane['S_h'] = S_h
    airplane['b_h'] = b_h
    airplane['xr_h'] = xr_h
    airplane['cr_h'] = cr_h
    airplane['xt_h'] = xt_h
    airplane['yt_h'] = yt_h
    airplane['zt_h'] = zt_h
    airplane['ct_h'] = ct_h
    airplane['xm_h'] = xm_h
    airplane['ym_h'] = ym_h
    airplane['zm_h'] = zm_h
    airplane['cm_h'] = cm_h
    airplane['S_v'] = S_v
    airplane['b_v'] = b_v
    airplane['xr_v'] = xr_v
    airplane['cr_v'] = cr_v
    airplane['xt_v'] = xt_v
    airplane['zt_v'] = zt_v
    airplane['ct_v'] = ct_v
    airplane['xm_v'] = xm_v
    airplane['zm_v'] = zm_v
    airplane['cm_v'] = cm_v

    # All variables are stored in the dictionary.
    # There is no need to return anything
    return None

#----------------------------------------

def aerodynamics(Mach, altitude, n_engines_failed, flap_def, slat_def,
                 lg_down, h_ground, W0_guess, airplane, method=2):
    '''
    method: 1 or 2 -> Method 1 applies a single friction coefficient
                      to the entire wetted area of the aircraft (based on Howe).
                      Method 2 is more refined since it computes friction and
                      form factors for each component.
    '''
    # c_flap_c_wing: extended total chord/ retracted total chord

    # Wetted areas from Torenbeek's Appendix B

    # Unpacking dictionary
    S_w = airplane['S_w']
    AR_w = airplane['AR_w']
    cr_w = airplane['cr_w']
    ct_w = airplane['ct_w']
    taper_w = airplane['taper_w']
    sweep_w = airplane['sweep_w']
    tcr_w = airplane['tcr_w']
    tct_w = airplane['tct_w']
    b_w = airplane['b_w']
    cm_w = airplane['cm_w']
    clmax_w = airplane['clmax_w']
    S_h = airplane['S_h']
    cr_h = airplane['cr_h']
    ct_h = airplane['ct_h']
    taper_h = airplane['taper_h']
    sweep_h = airplane['sweep_h']
    tcr_h = airplane['tcr_h']
    tct_h = airplane['tct_h']
    b_h = airplane['b_h']
    cm_h = airplane['cm_h']
    S_v = airplane['S_v']
    cr_v = airplane['cr_v']
    ct_v = airplane['ct_v']
    taper_v = airplane['taper_v']
    sweep_v = airplane['sweep_v']
    tcr_v = airplane['tcr_v']
    tct_v = airplane['tct_v']
    b_v = airplane['b_v']
    cm_v = airplane['cm_v']
    L_f = airplane['L_f']
    D_f = airplane['D_f']
    L_n = airplane['L_n']
    D_n = airplane['D_n']
    n_engines = airplane['n_engines']
    n_engines_under_wing = airplane['n_engines_under_wing']
    max_flap_def = max(airplane['TO_flap_def'], airplane['LD_flap_def'])
    flap_type = airplane['flap_type']
    c_flap_c_wing = airplane['c_flap_c_wing']
    b_flap_b_wing = airplane['b_flap_b_wing']
    max_slat_def = max(airplane['TO_slat_def'], airplane['LD_slat_def'])
    slat_type = airplane['slat_type']
    c_slat_c_wing = airplane['c_slat_c_wing']
    b_slat_b_wing = airplane['b_slat_b_wing']
    k_exc_drag = airplane['k_exc_drag']
    
    # Default rugosity value (smmoth paint from Raymer Tab 12.5)
    rugosity = 0.634e-5

    ### WING

    # Average t/c
    tc_avg = 0.5*(tcr_w + tct_w)

    #Exposed Area
    Sexp = S_w - cr_w*D_f

    #Wetted Area
    tau = tcr_w/tct_w
    Swet_w = 2*Sexp*(1 + 0.25*tcr_w*(1 + tau*taper_w)/(1 + taper_w))

    # Friction coefficient
    Cf_w = Cf_calc(Mach, altitude,
                   length = cm_w,
                   rugosity = rugosity,
                   k_lam = 0.05)

    # Form factor
    FF_w = FF_surface(Mach, tcr_w, tct_w, sweep_w, b_w, cr_w, ct_w, cm_w)

    # Interference factor
    Q_w = 1.0

    # Drag coefficient
    CD0_w = Cf_w*FF_w*Q_w*Swet_w/S_w

    ### HORIZONTAL TAIL

    #Exposed Area
    Sexp = S_h

    #Wetted Area
    tau = tcr_h/tct_h
    Swet_h = 2*Sexp*(1 + 0.25*tcr_h*(1 + tau*taper_h)/(1 + taper_h))

    # Friction coefficient
    Cf_h = Cf_calc(Mach, altitude,
                   length = cm_h,
                   rugosity = rugosity,
                   k_lam = 0.05)

    # Form factor
    FF_h = FF_surface(Mach, tcr_h, tct_h, sweep_h, b_h, cr_h, ct_h, cm_h)

    # Interference factor
    Q_h = 1.0

    # Drag coefficient
    CD0_h = Cf_h*FF_h*Q_h*Swet_h/S_w

    ### VERTICAL TAIL

    #Exposed Area
    Sexp = S_v

    #Wetted Area
    tau = tcr_v/tct_v
    Swet_v = 2*Sexp*(1 + 0.25*tcr_v*(1 + tau*taper_v)/(1 + taper_v))

    # Friction coefficient
    Cf_v = Cf_calc(Mach, altitude,
                   length = cm_v,
                   rugosity = rugosity,
                   k_lam = 0.05)

    # Form factor
    FF_v = FF_surface(Mach, tcr_v, tct_v, sweep_v, 2*b_v, cr_v, ct_v, cm_v)

    # Interference factor
    Q_v = 1.0

    # Drag coefficient
    CD0_v = Cf_v*FF_v*Q_v*Swet_v/S_w

    ### FUSELAGE

    # Wetted area
    lambda_fus = L_f/D_f
    Swet_f = np.pi*D_f*L_f*(1 - 2/lambda_fus)**(2.0/3.0)*(1 + 1/lambda_fus**2)

    # Friction coefficient
    Cf_f = Cf_calc(Mach, altitude,
                   length = L_f,
                   rugosity = rugosity,
                   k_lam = 0.05)

    # Form factor
    FF_f = 1 + 60/lambda_fus**3 + lambda_fus/400

    # Interference factor
    Q_f = 1.0

    # Drag coefficient
    CD0_f = Cf_f*FF_f*Q_f*Swet_f/S_w

    ### NACELLE

    # Wetted area (where we take the number of nacelles into account)
    Swet_n = n_engines*np.pi*D_n*L_n

    # Friction coefficient
    Cf_n = Cf_calc(Mach, altitude,
                   length = L_n,
                   rugosity = rugosity,
                   k_lam = 0.05)

    # Form factor
    lambda_n = L_n/D_n
    FF_n = 1 + 0.35/lambda_n

    # Interference factor
    Q_n = 1.2

    # Drag coefficient
    CD0_n = Cf_n*FF_n*Q_n*Swet_n/S_w

    ### VISCOUS DRAG

    if method == 1:
        
        # Total wetted area
        Swet = Swet_w + Swet_h + Swet_v + Swet_f + Swet_n
    
        # Wetted area ratio
        Sr = Swet/S_w
    
        # t/c correction
        tau = (Sr-2)/Sr + 1.9/Sr*(1 + 0.526*(4*tc_avg)**3)
    
        # Other parameters for jet aircraft
        Af = 0.93
        clam = 0.05
        Tf = 1.1
    
        # Friction coefficient (Howe Eq 6.13)
        Cfe = 0.005*(1-2*clam/Sr)*tau*(1 - 0.2*Mach + 0.12*(Mach*np.sqrt(np.cos(sweep_w))/(Af - tc_avg))**20)*Tf*S_w**(-0.1)
    
        # Viscous drag
        CD0 = Cfe*Swet/S_w

    elif method == 2:

        # Add all drag coefficients
        CD0 = CD0_w + CD0_h + CD0_v + CD0_f + CD0_n

    ### INDUCED

    # Oswald Factor (Howe Eq 6.14)
    f_taper = 0.005*(1 + 1.5*(taper_w - 0.6)**2)
    e = 1/(1 + 0.12*Mach**6)/(1 + (0.142 + AR_w*(10*tc_avg)**0.33*f_taper)/np.cos(sweep_w)**2 + 0.1*(3*n_engines_under_wing + 1)/(4 + AR_w)**0.8)

    # Induced drag term
    K = 1/np.pi/AR_w/e

    ### GROUND EFFECT
    if h_ground > 0:
        aux = 33*(h_ground/b_w)**1.5
        Kge = aux/(1+aux) # Raymer Eq. 12.61
        K = K*Kge

    ### Clean wing CLmax (Raymer Eq. 5.7)
    CLmax_clean = 0.9*clmax_w*np.cos(sweep_w)

    ### Flaps deflection
    ct_w = cr_w*taper_w
    if max_flap_def > 0.0:
        CD0_flap = 0.0023*b_flap_b_wing*flap_def*180/np.pi # Raymer Eq 12.37
        sweep_flap=geo_change_sweep(0.25, 2-c_flap_c_wing, sweep_w, b_w/2, cr_w, ct_w)
        if flap_type == 'plain':
            dclmax = 0.9
        elif flap_type == 'slotted':
            dclmax = 1.3
        elif flap_type == 'fowler':
            dclmax = 1.3*c_flap_c_wing
        elif flap_type == 'double slotted':
            dclmax = 1.6*c_flap_c_wing
        elif flap_type == 'triple slotted':
            dclmax = 1.9*c_flap_c_wing
        deltaCLmax_flap = dclmax*b_flap_b_wing*np.cos(sweep_flap)*flap_def/max_flap_def # Raymer Eq 12.21
    else:
        CD0_flap = 0.0
        deltaCLmax_flap = 0.0

    ### Slats deflection
    if max_slat_def > 0.0:
        CD0_slat = 0.0023*b_slat_b_wing*slat_def*180/np.pi # Raymer Eq 12.37
        sweep_slat=geo_change_sweep(0.25, c_slat_c_wing-1, sweep_w, b_w/2, cr_w, ct_w)
        if slat_type == 'fixed':
            dclmax = 0.2
        elif slat_type == 'flap':
            dclmax = 0.3
        elif slat_type == 'kruger':
            dclmax = 0.3
        elif slat_type == 'slat':
            dclmax = 0.4*c_slat_c_wing
        deltaCLmax_slat = dclmax*b_slat_b_wing*np.cos(sweep_slat)*slat_def/max_slat_def # Raymer Eq 12.21
    else:
        CD0_slat = 0.0
        deltaCLmax_slat = 0.0

    # Maximum lift
    CLmax = CLmax_clean + deltaCLmax_flap + deltaCLmax_slat

    ### Landing gear (ESDU)
    lg_factor = (0.57 - 0.26*flap_def/max_flap_def)*1e-3
    CD0_lg = lg_down*lg_factor*(W0_guess/gravity)**0.785/S_w

    ### Windmill engine
    #Vn_V = 0.42
    #CDwdm = (0.0785*D_n**2 + 1/(1 + 0.16*Mach**2)*np.pi/2*D_n**2*Vn_V*(1-Vn_V))/S_w
    #CD0_wdm = n_engines_failed*CDwdm
    CD0_wdm = n_engines_failed*0.3*np.pi/4*D_n**2/S_w # Raymer Eq 12.41

    # Add all drag values found so far
    CD0 = CD0 + CD0_flap + CD0_slat + CD0_lg + CD0_wdm

    ### Excrescence
    CD0_exc = CD0*k_exc_drag/(1-k_exc_drag)
    CD0 = CD0 + CD0_exc

    ### WAVE DRAG (Korn Equation)

    # The CL may get unreasonably high values for low Mach numbers.
    # This decreases Mdd, resulting in an unrealistic wave drag for
    # low Mach numbers. So we add another condition to discard any
    # wave drag below Mach 0.4

    if Mach > 0.4:

        # Estimate flight CL
        T,p,rho,mi = atmosphere(altitude)
        a = np.sqrt(1.4*287*T)
        V = a*Mach
        CL = 2*W0_guess/rho/V**2/S_w

        Mach_dd = 0.95/np.cos(sweep_w) - tc_avg/np.cos(sweep_w)**2 - CL/10/np.cos(sweep_w)**3
        Mach_crit = Mach_dd - (0.1/80)**(1/3)

        if (Mach > Mach_crit):
            CDw = 20*(Mach - Mach_crit)**4
        else:
            CDw = 0.0

    else:
        CDw = 0.0

    CD0 = CD0 + CDw

    # Update dictionary
    airplane['Swet_f'] = Swet_f

    return CD0, K, CLmax

#----------------------------------------

def engineTSFC(Mach, altitude, airplane):

    # Unpack dictionary
    BPR = airplane['BPR']
    Cbase = airplane['Cbase']

    ### ADD CODE FROM SECTION 3.3 HERE ###


    return C

#----------------------------------------

def empty_weight(W0_guess, T0_guess, airplane):

    # Unpack dictionary
    S_w = airplane['S_w']
    AR_w = airplane['AR_w']
    taper_w = airplane['taper_w']
    sweep_w = airplane['sweep_w']
    xm_w = airplane['xm_w']
    cm_w = airplane['cm_w']
    tcr_w = airplane['tcr_w']
    S_h = airplane['S_h']
    xm_h = airplane['xm_h']
    cm_h = airplane['cm_h']
    S_v = airplane['S_v']
    xm_v = airplane['xm_v']
    cm_v = airplane['cm_v']
    L_f = airplane['L_f']
    Swet_f = airplane['Swet_f']
    n_engines = airplane['n_engines']
    BPR = airplane['BPR']
    x_n = airplane['x_n']
    L_n = airplane['L_n']
    x_nlg = airplane['x_nlg']
    x_mlg = airplane['x_mlg']

    ### ADD CODE FROM SECTION 3.4 HERE ###


    # Update dictionary
    airplane['W_w'] = W_w
    airplane['W_h'] = W_h
    airplane['W_v'] = W_v
    airplane['W_f'] = W_f
    airplane['W_nlg'] = W_nlg
    airplane['W_mlg'] = W_mlg
    airplane['W_eng'] = W_eng_installed
    airplane['W_allelse'] = W_allelse

    return We, xcg_e

#----------------------------------------

def fuel_weight(W0_guess, airplane):

    # Unpacking dictionary
    S_w = airplane['S_w']
    altitude_cruise = airplane['altitude_cruise']
    Mach_cruise = airplane['Mach_cruise']
    range_cruise = airplane['range_cruise']
    loiter_time = airplane['loiter_time']
    altitude_altcruise = airplane['altitude_altcruise']
    Mach_altcruise = airplane['Mach_altcruise']
    range_altcruise = airplane['range_altcruise']

    ### ADD CODE FROM SECTION 3.5 HERE ###


    return Wf, Mf_cruise

#----------------------------------------

def weight(W0_guess, T0_guess, airplane):

    # Unpacking dictionary
    W_payload = airplane['W_payload']
    W_crew = airplane['W_crew']

    # Set iterator
    delta = 1000

    while abs(delta) > 10:
        pass
        ### ADD CODE FROM SECTION 3.6.4 HERE ###


    return W0, We, Wf, Mf_cruise, xcg_e

#----------------------------------------

def performance(W0, Mf_cruise, airplane):

    '''
    This function computes the required thrust and wing areas
    required to meet takeoff, landing, climb, and cruise requirements.

    OUTPUTS:
    T0: real -> Total thrust required to meet all mission phases
    S_wlan: real -> Wing area required for landing. The wing area (S_w) should
                    be greater than this value.
    '''

    # Unpacking dictionary
    S_w = airplane['S_w']
    n_engines = airplane['n_engines']
    BPR = airplane['BPR']
    TO_flap_def = airplane['TO_flap_def']
    LD_flap_def = airplane['LD_flap_def']
    TO_slat_def = airplane['TO_slat_def']
    LD_slat_def = airplane['LD_slat_def']
    h_ground = airplane['h_ground']
    altitude_takeoff = airplane['altitude_takeoff']
    distance_takeoff = airplane['distance_takeoff']
    altitude_landing = airplane['altitude_landing']
    distance_landing = airplane['distance_landing']
    MLW_frac = airplane['MLW_frac']
    altitude_cruise = airplane['altitude_cruise']
    Mach_cruise = airplane['Mach_cruise']
    
    ### ADD CODE FROM SECTION 3.7.3 TO SECTION 3.7.5 HERE ###


    ### CLIMB

    # Define standard function for climb analysis
    def climb_analysis(grad, Ks, altitude, CLmax_guess,
                       lg_down, h_ground_climb, flap_def, slat_def, n_engines_failed, Mf,
                       kT):

        '''
        We need a guess for CLmax just to get an approximate drag polar for
        speed computation. We will get the correct CLmax from the aerodynamics module

        kT: Thrust decay factor (e.g. use 0.94 for maximum continuous thrust)
        '''

        ### ADD CODE FROM SECTION 3.7.6 HERE ###


        return T0

    ### CONTINUE THE CODE FROM SECTION 3.7.6 HERE ###

    ### ADD CODE FROM SECTION 3.7.7 HERE ###

    return T0, T0vec, deltaS_wlan, CLmaxTO

#----------------------------------------

def thrust_matching(W0_guess, T0_guess, airplane):

    ### ADD CODE FROM SECTION 3.8 HERE ###

    # Update dictionary
    airplane['W0'] = W0
    airplane['We'] = We
    airplane['Wf'] = Wf
    airplane['xcg_e'] = xcg_e
    airplane['T0'] = T0
    airplane['T0vec'] = T0vec
    airplane['deltaS_wlan'] = deltaS_wlan
    airplane['CLmaxTO'] = CLmaxTO

    # Return
    return None

#----------------------------------------

def balance(airplane):

    # Unpack dictionary
    W0 = airplane['W0']
    W_payload = airplane['W_payload']
    xcg_payload = airplane['xcg_payload']
    W_crew = airplane['W_crew']
    xcg_crew = airplane['xcg_crew']
    We = airplane['We']
    xcg_e = airplane['xcg_e']
    Wf = airplane['Wf']
    Mach_cruise = airplane['Mach_cruise']
    S_w = airplane['S_w']
    AR_w = airplane['AR_w']
    sweep_w = airplane['sweep_w']
    b_w = airplane['b_w']
    xr_w = airplane['xr_w']
    cr_w = airplane['cr_w']
    ct_w = airplane['ct_w']
    xm_w = airplane['xm_w']
    cm_w = airplane['cm_w']
    tcr_w = airplane['tcr_w']
    tct_w = airplane['tct_w']
    c_tank_c_w = airplane['c_tank_c_w']
    x_tank_c_w = airplane['x_tank_c_w']
    S_h = airplane['S_h']
    AR_h = airplane['AR_h']
    sweep_h = airplane['sweep_h']
    b_h = airplane['b_h']
    cr_h = airplane['cr_h']
    ct_h = airplane['ct_h']
    xm_h = airplane['xm_h']
    cm_h = airplane['cm_h']
    eta_h = airplane['eta_h']
    Cvt = airplane['Cvt']
    L_f = airplane['L_f']
    D_f = airplane['D_f']
    y_n = airplane['y_n']
    T0 = airplane['T0']
    n_engines = airplane['n_engines']
    CLmaxTO = airplane['CLmaxTO']
    rho_f = airplane['rho_f']
    
    ### ADD CODE FROM SECTION 3.9 HERE ###

    # Update dictionary
    airplane['xcg_fwd'] = xcg_fwd
    airplane['xcg_aft'] = xcg_aft
    airplane['xnp'] = xnp
    airplane['SM_fwd'] = SM_fwd
    airplane['SM_aft'] = SM_aft
    airplane['b_tank_b_w'] = b_tank_b_w
    airplane['CLv'] = CLv

    return None

#----------------------------------------

def landing_gear(airplane):

    # Unpack dictionary
    x_nlg = airplane['x_nlg']
    x_mlg = airplane['x_mlg']
    y_mlg = airplane['y_mlg']
    z_lg = airplane['z_lg']
    xcg_fwd = airplane['xcg_fwd']
    xcg_aft = airplane['xcg_aft']
    x_tailstrike = airplane['x_tailstrike']
    z_tailstrike = airplane['z_tailstrike']

    ### ADD CODE FROM SECTION 3.10 HERE ###

    # Update dictionary
    airplane['frac_nlg_fwd'] = frac_nlg_fwd
    airplane['frac_nlg_aft'] = frac_nlg_aft
    airplane['alpha_tipback'] = alpha_tipback
    airplane['alpha_tailstrike'] = alpha_tailstrike
    airplane['phi_overturn'] = phi_overturn

    return None

#========================================
# AUXILIARY FUNCTIONS

def plot3d(airplane, figname='3dview.png'):

    '''
    This function generates a 3D plot of the aircraft
    '''

    # Unpack dictionary
    xr_w = airplane['xr_w']
    zr_w = airplane['zr_w']
    cr_w = airplane['cr_w']
    xt_w = airplane['xt_w']
    yt_w = airplane['yt_w']
    zt_w = airplane['zt_w']
    ct_w = airplane['ct_w']
    xm_w = airplane['xm_w']
    ym_w = airplane['ym_w']
    zm_w = airplane['zm_w']
    cm_w = airplane['cm_w']
    xr_h = airplane['xr_h']
    zr_h = airplane['zr_h']
    cr_h = airplane['cr_h']
    xt_h = airplane['xt_h']
    yt_h = airplane['yt_h']
    zt_h = airplane['zt_h']
    ct_h = airplane['ct_h']
    xm_h = airplane['xm_h']
    ym_h = airplane['ym_h']
    zm_h = airplane['zm_h']
    cm_h = airplane['cm_h']
    xr_v = airplane['xr_v']
    zr_v = airplane['zr_v']
    cr_v = airplane['cr_v']
    xt_v = airplane['xt_v']
    zt_v = airplane['zt_v']
    ct_v = airplane['ct_v']
    xm_v = airplane['xm_v']
    zm_v = airplane['zm_v']
    cm_v = airplane['cm_v']
    L_f = airplane['L_f']
    D_f = airplane['D_f']
    x_n = airplane['x_n']
    y_n = airplane['y_n']
    z_n = airplane['z_n']
    L_n = airplane['L_n']
    D_n = airplane['D_n']

    # Optional arguments

    if 'xcg_fwd' in airplane:
        xcg_fwd = airplane['xcg_fwd']
        xcg_aft = airplane['xcg_aft']
    else:
        xcg_fwd = None
        xcg_aft = None

    if 'xnp' in airplane:
        xnp = airplane['xnp']
    else:
        xnp = None
           
    ### PLOT

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    # Wing plot
    ax.plot([xr_w, xt_w, xt_w+ct_w, xr_w+cr_w, xt_w+ct_w, xt_w, xr_w],
            [0.0, yt_w, yt_w, 0.0, -yt_w, -yt_w, 0.0],
            [zr_w, zt_w, zt_w, zr_w, zt_w, zt_w, zr_w])

    # HT plot
    ax.plot([xr_h, xt_h, xt_h+ct_h, xr_h+cr_h, xt_h+ct_h, xt_h, xr_h],
            [0.0, yt_h, yt_h, 0.0, -yt_h, -yt_h, 0.0],
            [zr_h, zt_h, zt_h, zr_h, zt_h, zt_h, zr_h])

    # VT plot
    ax.plot([xr_v, xt_v, xt_v+ct_v, xr_v+cr_v, xr_v],
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [zr_v, zt_v, zt_v, zr_v, zr_v])

    # Fuselage plot
    ax.plot([0.0, L_f],
            [0.0, 0.0],
            [0.0, 0.0])

    # Nacelle 1 plot
    ax.plot([x_n, x_n+L_n],
            [y_n, y_n],
            [z_n, z_n])

    # Nacelle 2 plot
    ax.plot([x_n, x_n+L_n],
            [-y_n, -y_n],
            [z_n, z_n])

    # Mean aerodynamic chords
    ax.plot([xm_w, xm_w+cm_w],
            [ym_w, ym_w],
            [zm_w, zm_w],'g')
    ax.plot([xm_h, xm_h+cm_h],
            [ym_h, ym_h],
            [zm_h, zm_h],'g')
    ax.plot([xm_v, xm_v+cm_v],
            [0.0, 0.0],
            [zm_v, zm_v],'g')

    # Center of gravity
    if xcg_fwd is not None:
        ax.plot([xcg_fwd, xcg_aft],
                [0.0, 0.0],
                [0.0, 0.0],'o')

    # Neutral point
    if xnp is not None:
        ax.plot([xnp, xnp],
                [0.0, 0.0],
                [0.0, 0.0],'o')

    # Create cubic bounding box to simulate equal aspect ratio
    X = np.array([xr_w, xt_h+ct_h, xt_v+ct_v])
    Y = np.array([-yt_w, yt_w])
    Z = np.array([zt_w, zt_h, zt_v])
    max_range = np.array([X.max()-X.min(), Y.max()-Y.min(), Z.max()-Z.min()]).max()
    Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(X.max()+X.min())
    Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(Y.max()+Y.min())
    Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(Z.max()+Z.min())

    # Comment or uncomment following both lines to test the fake bounding box:
    for xb, yb, zb in zip(Xb, Yb, Zb):
        ax.plot([xb], [yb], [zb], 'w')

    # Set initial point of view
    ax.view_init(45,-135)

    # Save figure
    plt.show()
    fig.savefig(figname,dpi=300)

#----------------------------------------

def atmosphere(z, Tba=288.15):

    '''
    Funçao que retorna a Temperatura, Pressao e Densidade para uma determinada
    altitude z [m]. Essa funçao usa o modelo padrao de atmosfera para a
    temperatura no solo de Tba.
    '''

    # Zbase (so para referencia)
    # 0 11019.1 20063.1 32161.9 47350.1 50396.4

    # DEFINING CONSTANTS
    # Earth radius
    r = 6356766
    # gravity
    g0 = 9.80665
    # air gas constant
    R = 287.05287
    # layer boundaries
    Ht = [0, 11000, 20000, 32000, 47000, 50000]
    # temperature slope in each layer
    A = [-6.5e-3, 0, 1e-3, 2.8e-3, 0]
    # pressure at the base of each layer
    pb = [101325, 22632, 5474.87, 868.014, 110.906]
    # temperature at the base of each layer
    Tstdb = [288.15, 216.65, 216.65, 228.65, 270.65];
    # temperature correction
    Tb = Tba-Tstdb[0]
    # air viscosity
    mi0 = 18.27e-6 # [Pa s]
    T0 = 291.15 # [K]
    C = 120 # [K]

    # geopotential altitude
    H = r*z/(r+z)

    # selecting layer
    if H < Ht[0]:
        raise ValueError('Under sealevel')
    elif H <= Ht[1]:
        i = 0
    elif H <= Ht[2]:
        i = 1
    elif H <= Ht[3]:
        i = 2
    elif H <= Ht[4]:
        i = 3
    elif H <= Ht[5]:
        i = 4
    else:
        raise ValueError('Altitude beyond model boundaries')

    # Calculating temperature
    T = Tstdb[i]+A[i]*(H-Ht[i])+Tb

    # Calculating pressure
    if A[i] == 0:
        p = pb[i]*np.exp(-g0*(H-Ht[i])/R/(Tstdb[i]+Tb))
    else:
        p = pb[i]*(T/(Tstdb[i]+Tb))**(-g0/A[i]/R)

    # Calculating density
    rho = p/R/T

    # Calculating viscosity with Sutherland's Formula
    mi=mi0*(T0+C)/(T+C)*(T/T0)**(1.5)

    return T,p,rho,mi

#----------------------------------------

def geo_change_sweep(x,y,sweep_x,panel_length,chord_root,chord_tip):

    '''
    This function converts sweep computed at chord fraction x into
    sweep measured at chord fraction y
    (x and y should be between 0 (leading edge) and 1 (trailing edge).
    '''

    sweep_y=sweep_x+np.arctan((x-y)*(chord_root-chord_tip)/panel_length)

    return sweep_y

#----------------------------------------

def Cf_calc(Mach, altitude, length, rugosity, k_lam, Tba=288.15):
    '''
    This function computes the flat plate friction coefficient
    for a given Reynolds number while taking transition into account

    k_lam: float -> Fraction of the length (from 0 to 1) where
                    transition occurs
    '''
    
    # Dados atmosféricos
    T, p, rho, mi = atmosphere(altitude, Tba)

    # Velocidade
    v = np.sqrt(1.4*287*T)*Mach

    # Reynolds na transição
    Re_conv = rho*v*k_lam*length/mi
    Re_rug = 38.21*(k_lam*length/rugosity)**1.053
    Re_trans = min(Re_conv, Re_rug)

    # Reynolds no fim
    Re_conv = rho*v*length/mi
    Re_rug = 38.21*(length/rugosity)**1.053
    Re_fim = min(Re_conv, Re_rug)

    # Coeficientes de fricção
    # Laminar na transição
    Cf1 = 1.328/np.sqrt(Re_trans)

    # Turbulento na transição
    Cf2 = 0.455/(np.log10(Re_trans)**2.58*(1+0.144*Mach**2)**0.65)

    # Turbulento no fim
    Cf3 = 0.455/(np.log10(Re_fim)**2.58*(1+0.144*Mach**2)**0.65)

    # Média
    Cf = (Cf1 - Cf2)*k_lam + Cf3

    return Cf

#----------------------------------------

def FF_surface(Mach, tcr, tct, sweep, b, cr, ct, cm, x_c_max_tc=0.4):
    '''
    This function computes the form factor for lifting surfaces

    INPUTS

    tcr: float -> Thickness/chord ratio at the root
    tct: float -> Thickness/chord ratio at the tip
    sweep: float -> Quarter-chord sweep angle [rad]
    b: float -> Wing span (considering both sides. Double this value for vertical tails if necessary)
    cr: float -> Root chord
    ct: float -> Tip chord
    cm: float -> Mean aerodynamic chord
    x_c_max_tc: float -> Chord fraction with maximum thickness
    '''

    # Average chord fraction
    t_c = (tcr + tct)/2

    # Swee at maximum thickness position
    sweep_maxtc=geo_change_sweep(0.25, x_c_max_tc, sweep, b/2, cr, ct)

    # Form factor
    FF = 1.34*Mach**0.18*np.cos(sweep_maxtc)**0.28*(1 + 0.6*t_c/x_c_max_tc/cm + 100*(t_c)**4)

    return FF

#----------------------------------------

def standard_airplane():
    '''
    The standard parameters refer to the Fokker 100, but they could be redefined for
    any new aircraft.
    '''

    airplane = {'S_w' : 93.5, # Wing area [m2]
                'AR_w' : 8.43,  # Wing aspect ratio
                'taper_w' : 0.235, # Wing taper ratio
                'sweep_w' : 17.45*np.pi/180, # Wing sweep [rad]
                'dihedral_w' : 5*np.pi/180, # Wing dihedral [rad]
                'xr_w' : 13.5, # Longitudinal position of the wing (with respect to the fuselage nose) [m]
                'zr_w' : 0.0, # Vertical position of the wing (with respect to the fuselage nose) [m]
                'tcr_w' : 0.123, # t/c of the root section of the wing
                'tct_w' : 0.096, # t/c of the tip section of the wing
                
                'Cht' : 0.94, # Horizontal tail volume coefficient
                'Lc_h' : 4.83, # Non-dimensional lever of the horizontal tail (lever/wing_mac)
                'AR_h' : 4.64, # HT aspect ratio
                'taper_h' : 0.39, # HT taper ratio
                'sweep_h' : 26*np.pi/180, # HT sweep [rad]
                'dihedral_h' : 2*np.pi/180, # HT dihedral [rad]
                'zr_h' : 4.359, # Vertical position of the HT [m]
                'tcr_h' : 0.1, # t/c of the root section of the HT
                'tct_h' : 0.1, # t/c of the tip section of the HT
                'eta_h' : 1.0, # Dynamic pressure factor of the HT
                
                'Cvt' : 0.088, # Vertical tail volume coefficient
                'Lb_v' : 0.55, # Non-dimensional lever of the vertical tail (lever/wing_span)
                'AR_v' : 1.27, # VT aspect ratio
                'taper_v' : 0.74, # VT taper ratio
                'sweep_v' : 41*np.pi/180, # VT sweep [rad]
                'zr_v' : 0.0, # Vertical position of the VT [m]
                'tcr_v' : 0.1, # t/c of the root section of the VT
                'tct_v' : 0.1, # t/c of the tip section of the VT
                
                'L_f' : 32.5, # Fuselage length [m]
                'D_f' : 3.3, # Fuselage diameter [m]
                
                'x_n' : 23.2, # Longitudinal position of the nacelle frontal face [m]
                'y_n' : 2.6, # Lateral position of the nacelle centerline [m]
                'z_n' : 0.0, # Vertical position of the nacelle centerline [m]
                'L_n' : 4.3, # Nacelle length [m]
                'D_n' : 1.5, # Nacelle diameter [m]
                
                'n_engines' : 2, # Number of engines
                'n_engines_under_wing' : 0, # Number of engines installed under the wing
                'BPR' : 3.04, # Engine bypass ratio
                'Cbase' : None, # Base engine TSFC [1/s] (use 'None' for Howe's values)
                
                'x_nlg' : 3.6, # Longitudinal position of the nose landing gear [m]
                'x_mlg' : 17.8, # Longitudinal position of the main landing gear [m]
                'y_mlg' : 2.47, # Lateral position of the main landing gear [m]
                'z_lg' : -2.0, # Vertical position of the landing gear [m]
                'x_tailstrike' : 23.68, # Longitudinal position of critical tailstrike point [m]
                'z_tailstrike' : -0.84, # Vertical position of critical tailstrike point [m]
                
                'c_tank_c_w' : 0.4, # Fraction of the wing chord occupied by the fuel tank
                'x_tank_c_w' : 0.2, # Fraction of the wing chord where fuel tank starts
                
                'clmax_w' : 2.1, # Maximum lift coefficient of wing airfoil
    
                'TO_flap_def' : 20*np.pi/180, # Takeoff flap deflection [rad]
                'LD_flap_def' : 40*np.pi/180, # Landing flap deflection [rad]
                'flap_type' : 'double slotted',  # Flap type
                'c_flap_c_wing' : 1.2, # chord_with_deflected_flaps/chord_with_retracted_flaps
                'b_flap_b_wing' : 0.6, # Fraction of the wing span occupied by flaps
                
                'TO_slat_def' : 0*np.pi/180, # Takeoff slat deflection [rad]
                'LD_slat_def' : 0*np.pi/180, # Landing slat deflection [rad]
                'slat_type' : 'slat', # Slat type
                'c_slat_c_wing' : 1.05, # chord_with_deflected_slats/chord_with_retracted_slats
                'b_slat_b_wing' : 0.75, # Fraction of the wing span occupied by slats
                
                'h_ground' : 35.0*ft2m, # Distance to the ground for ground effect computation [m]
                'k_exc_drag' : 0.03, # Excrescence drag factor
                
                'altitude_takeoff' : 0.0, # Altitude for takeoff computation [m]
                'distance_takeoff' : 1800.0, # Required takeoff distance [m]
                
                'altitude_landing' : 0.0, # Altitude for landing computation [m]
                'distance_landing' : 1800.0, # Required landing distance [m] (The actual Fokker100 distance is 1350 m but it is very restrictive compared to the historical regression. Therefore I kept the same TO distance since the aircraft should takeoff and land at the same runway)
                'MLW_frac' : 38300/41500, # Max Landing Weight / Max Takeoff Weight
                
                'altitude_cruise' : 35000*ft2m, # Cruise altitude [m]
                'Mach_cruise' : 0.73, # Cruise Mach number
                'range_cruise' : 1200*nm2m, # Cruise range [m]
                
                'loiter_time' : 45*60, # Loiter time [s]
                
                'altitude_altcruise' : 4572, # Alternative cruise altitude [m]
                'Mach_altcruise' : 0.4, # Alternative cruise Mach number
                'range_altcruise' : 200*nm2m, # Alternative cruise range [m]
                
                'W_payload' : 107*91*gravity, # Payload weight [N]
                'xcg_payload' : 14.4, # Longitudinal position of the Payload center of gravity [m]
                
                'W_crew' : 5*91*gravity, # Crew weight [N]
                'xcg_crew' : 2.5, # Longitudinal position of the Crew center of gravity [m]
                
                'rho_f' : 804, # Fuel density kg/m3 (This is Jet A-1)
                }

    return airplane

