import designTool as dt
import numpy as np
import matplotlib.pyplot as plt

skysonic = {'AR_h': 4.19,
            'AR_v': 1.22,
            'AR_w': 6.92,
            'BPR': 3.04,
            'Cbase': None,
            'Cht': 0.8,
            'Cvt': 0.081,
            'D_f': 2.79,
            'D_n': 1.5,
            'LD_flap_def': 0.6981317007977318,
            'LD_slat_def': 0.0,
            'L_f': 23.19,
            'L_n': 4.3,
            'Lb_v': 0.448,
            'Lc_h': 3.5,  # calculado e ajustado com base no ERJ-145 e CRJ-200
            'MLW_frac': 0.9228915662650602,
            'Mach_altcruise': 0.4,  # não sabemos o que é
            'Mach_cruise': 0.75,
            'S_w': 68.89,               #!
            'Swet_f': 1,
            'TO_flap_def': 0.3490658503988659,
            'TO_slat_def': 0.0,
            'W_crew': 1240*dt.lb2N,
            'W_payload': 15036.67*dt.lb2N,
            'altitude_altcruise': 4572,
            'altitude_cruise': 12192,
            'altitude_landing': 0.0,
            'altitude_takeoff': 0.0,
            'b_flap_b_wing': 0.6,
            'b_slat_b_wing': 0.75,
            'c_flap_c_wing': 1.2,
            'c_slat_c_wing': 1.05,
            'c_tank_c_w': 0.515,          #! was 0.4, worked with 0.515
            'clmax_w': 2.1,
            #!
            'cm_h': 2.107457619636192,
            'cm_v': 3.4576757510555542,
            'cm_w': 3.756317488774531,
            'cr_h': 2.849393124273043,
            'cr_v': 3.944978890651773,
            'cr_w': 5.3933059334262,
            'ct_h': 1.1112633184664868,
            'ct_v': 2.919284379082312,
            'ct_w': 1.267426894355157,
            #!
            'dihedral_h': np.deg2rad(0),
            'dihedral_w': np.deg2rad(3.5),
            'distance_landing': 1800.0,
            'distance_takeoff': 1800.0,
            'eta_h': 1.0,
            'flap_type': 'double slotted',
            # distância do solo para que seja computado o efeito de solo
            'h_ground': 10.668000000000001,
            'k_exc_drag': 0.03,
            'loiter_time': 45*60 + ((dt.nm2m * 200) / (0.75 * 343)),
            'n_engines': 2,
            'n_engines_under_wing': 0,
            'range_altcruise': 370400.0,
            'range_cruise': dt.nm2m * 2000,
            'rho_f': 804,
            'slat_type': 'slat',
            'sweep_h': np.deg2rad(25),
            # alterado para encaixar melhor a EV e EH
            'sweep_v': np.deg2rad(40.5),
            'sweep_w': np.deg2rad(21),
            'taper_h': 0.52,
            'taper_v': 0.71,
            'taper_w': 0.3,
            'tcr_h': 0.1,
            'tcr_v': 0.1,
            'tcr_w': 0.123,
            'tct_h': 0.1,
            'tct_v': 0.1,
            'tct_w': 0.096,
            'x_mlg': 14.12,                        #! Testes   14.12
            'x_n': 16.76,
            'x_nlg': 2,                         #!
            'x_tailstrike': 23.68,              #! 
            'x_tank_c_w': 0.2,
            'xcg_crew': 2.5,
            'xcg_payload': 10,                  #! was 14.4, estimated 18.17
            'xr_w': 9.54,                       #! estimado com base em dados do CRJ-200
            'y_mlg': 2.47,
            'y_n': 2.6,
            'z_lg': -2.0,
            'z_n': 0.0,
            'z_tailstrike': -0.84,              #!
            'zr_h': 3.784,
            'zr_v': 0.0,
            'zr_w': 0.0}
skysonic['frac_xcg_ae'] = None
W0_guess = 467500
T0_guess = 140250

dt.analyze(airplane=skysonic, print_log=True, plot=True, W0_guess=W0_guess, T0_guess=T0_guess, review=True)
'''
frac_xcg_AllElses = np.linspace(0.4, 0.5, 100)
deltaS_wlan = np.zeros(len(frac_xcg_AllElses))
SM_aft = np.zeros(len(frac_xcg_AllElses))
SM_fwd = np.zeros(len(frac_xcg_AllElses))
CLv = np.zeros(len(frac_xcg_AllElses))
frac_nlg_fwd = np.zeros(len(frac_xcg_AllElses))
frac_nlg_aft = np.zeros(len(frac_xcg_AllElses))
alpha_tipback = np.zeros(len(frac_xcg_AllElses))
alpha_tailstrike = np.zeros(len(frac_xcg_AllElses))
phi_overturn = np.zeros(len(frac_xcg_AllElses))
T0_vecs = np.zeros((len(frac_xcg_AllElses), 8))

for i, fracs in enumerate(frac_xcg_AllElses):
    skysonic['frac_xcg_ae'] = fracs
    dt.analyze(skysonic, print_log=False, plot=False, W0_guess=W0_guess, T0_guess=T0_guess, review=False)  
    deltaS_wlan[i] = skysonic['deltaS_wlan']
    SM_aft[i] = skysonic['SM_aft']
    SM_fwd[i] = skysonic['SM_fwd']
    CLv[i] = skysonic['CLv']
    frac_nlg_fwd[i] = skysonic['frac_nlg_fwd']
    frac_nlg_aft[i] = skysonic['frac_nlg_aft']
    alpha_tipback[i] = skysonic['alpha_tipback']
    alpha_tailstrike[i] = skysonic['alpha_tailstrike']
    phi_overturn[i] = skysonic['phi_overturn']
    T0_vecs[i] = skysonic['T0vec'] 
    
plt.figure()
plt.plot(frac_xcg_AllElses*100, deltaS_wlan)
# plt.scatter(airplane['frac_xcg_ae']*100, airplane['deltaS_wlan'], marker='*', color='green', s=100)
plt.title('deltaS_wlan')
plt.axhline(0, linestyle='--', color = 'red', label='Limite')
plt.xlabel('Fração Xcg_ae [%]')
plt.legend()
    
plt.figure()
plt.plot(frac_xcg_AllElses*100, SM_aft)
# plt.scatter(airplane['frac_xcg_ae']*100, airplane['SM_aft'], marker='*', color='green', s=100)
plt.title('SM_aft')
plt.axhline(0.05, linestyle='--', color = 'red', label='Limite')
plt.xlabel('Fração Xcg_ae [%]')
plt.legend()

plt.figure()
plt.plot(frac_xcg_AllElses*100, SM_fwd)
plt.title('SM_fwd')
plt.axhline(0.30, linestyle='--', color = 'red', label='Limite')
# plt.scatter(airplane['frac_xcg_ae']*100, airplane['SM_fwd'], marker='*', color='green', s=100)
plt.xlabel('Fração Xcg_ae [%]')
plt.legend()

plt.figure()
plt.plot(frac_xcg_AllElses*100, CLv)
plt.title('CLv')
plt.axhline(0.75, linestyle='--', color = 'red', label='Limite')
# plt.scatter(airplane['frac_xcg_ae']*100, airplane['CLv'], marker='*', color='green', s=100)
plt.xlabel('Fração Xcg_ae [%]')
plt.legend()

plt.figure()
plt.plot(frac_xcg_AllElses*100, frac_nlg_fwd)
plt.title('frac_nlg_fwd')
plt.axhline(0.18, linestyle='--', color = 'red', label='Limite')
# plt.scatter(airplane['frac_xcg_ae']*100, airplane['frac_nlg_fwd'], marker='*', color='green', s=100)
plt.xlabel('Fração Xcg_ae [%]')
plt.legend()

plt.figure()
plt.plot(frac_xcg_AllElses*100, frac_nlg_aft)
plt.title('frac_nlg_aft [%]')
plt.axhline(0.05, linestyle='--', color = 'red', label='Limite')
# plt.scatter(airplane['frac_xcg_ae']*100, airplane['frac_nlg_aft'], marker='*', color='green', s=100)
plt.xlabel('Fração Xcg_ae [%]')
plt.legend()


plt.figure()
plt.plot(frac_xcg_AllElses*100, np.rad2deg(alpha_tipback))
# plt.scatter(airplane['frac_xcg_ae']*100, np.rad2deg(airplane['alpha_tipback']), marker='*', color='green', s=100)
plt.title('alpha_tipback')
plt.axhline(15, linestyle='--', color = 'red', label='Limite')
plt.xlabel('Fração Xcg_ae [%]')
plt.legend()


plt.figure()
plt.plot(frac_xcg_AllElses*100, np.rad2deg(alpha_tailstrike))
# plt.scatter(airplane['frac_xcg_ae']*100, np.rad2deg(airplane['alpha_tailstrike']), marker='*', color='green', s=100)
plt.title('alpha_tailstrike')
plt.axhline(10, linestyle='--', color = 'red', label='Limite')
plt.xlabel('Fração Xcg_ae [%]')
plt.legend()


plt.figure()
plt.plot(frac_xcg_AllElses*100, np.rad2deg(phi_overturn))
# plt.scatter(airplane['frac_xcg_ae']*100, np.rad2deg(airplane['phi_overturn']), marker='*', color='green', s=100)
plt.title('phi_overturn')
plt.axhline(63, linestyle='--', color = 'red', label='Limite')
plt.xlabel('Fração Xcg_ae [%]')
plt.legend()
plt.show()
'''