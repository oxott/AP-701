import designTool_david as dt
import numpy as np
import matplotlib.pyplot as plt
import pprint

# 3.9

airplane = { 'AR_h': 4.64,
'AR_v': 1.27,
'AR_w': 8.43,
'BPR': 3.04,
'CLmaxTO': 2.373034560798506,
'Cbase': None,
'Cht': 0.94,
'Cvt': 0.088,
'D_f': 3.3,
'D_n': 1.5,
'LD_flap_def': 0.6981317007977318,
'LD_slat_def': 0.0,
'L_f': 32.5,
'L_n': 4.3,
'Lb_v': 0.55,
'Lc_h': 4.83,
'MLW_frac': 0.9228915662650602,
'Mach_altcruise': 0.4,
'Mach_cruise': 0.73,
'S_h': 18.196687370600415,
'S_v': 14.959999999999999,
'S_w': 93.5,
'Swet_f': 292.60345689585,
'T0': 116905.51958773875,
'T0vec':   [103943.7814999834,
            85542.04301785328,
            91289.98209134683,
            96929.35159443578,
            108591.15496267364,
            74004.46014201488,
            62300.66945993366,
            111338.5900835607],
'TO_flap_def': 0.3490658503988659,
'TO_slat_def': 0.0,
'W0': 417029.8478522959,
'W_allelse': 70896.60884061972,
'W_crew': 4463.55,
'W_eng': 28712.35646954906,
'W_f': 68890.55789155893,
'W_h': 4819.756583850933,
'W_mlg': 15242.77090073324,
'W_nlg': 2689.9007471882187,
'W_payload': 95519.97,
'W_v': 3962.4552,
'W_w': 32692.945285195954,
'We': 227907.35191869608,
'Wf': 89138.97593359985,
'altitude_altcruise': 4572,
'altitude_cruise': 10668.0,
'altitude_landing': 0.0,
'altitude_takeoff': 0.0,
'b_flap_b_wing': 0.6,
'b_h': 9.18872294715571,
'b_slat_b_wing': 0.75,
'b_v': 4.358807176281144,
'b_w': 28.074988869098416,
'c_flap_c_wing': 1.2,
'c_slat_c_wing': 1.05,
'c_tank_c_w': 0.4,
'clmax_w': 2.1,
'cm_h': 2.107457619636192,
'cm_v': 3.4576757510555542,
'cm_w': 3.756317488774531,
'cr_h': 2.849393124273043,
'cr_v': 3.944978890651773,
'cr_w': 5.3933059334262,
'ct_h': 1.1112633184664868,
'ct_v': 2.919284379082312,
'ct_w': 1.267426894355157,
'deltaS_wlan': 24.218094249619313,
'dihedral_h': 0.03490658503988659,
'dihedral_w': 0.08726646259971647,
'distance_landing': 1800.0,
'distance_takeoff': 1800.0,
'eta_h': 1.0,
'flap_type': 'double slotted',
'h_ground': 10.668000000000001,
'k_exc_drag': 0.03,
'loiter_time': 2700,
'n_engines': 2,
'n_engines_under_wing': 0,
'range_altcruise': 370400.0,
'range_cruise': 2222400.0,
'rho_f': 804,
'slat_type': 'slat',
'sweep_h': 0.4537856055185257,
'sweep_v': 0.715584993317675,
'sweep_w': 0.3045599544730105,
'taper_h': 0.39,
'taper_v': 0.74,
'taper_w': 0.235,
'tcr_h': 0.1,
'tcr_v': 0.1,
'tcr_w': 0.123,
'tct_h': 0.1,
'tct_v': 0.1,
'tct_w': 0.096,
'x_mlg': 17.8,
'x_n': 23.2,
'x_nlg': 3.6,
'x_tailstrike': 23.68,
'x_tank_c_w': 0.2,
'xcg_crew': 2.5,
'xcg_e': 17.166310753037887,
'xcg_payload': 14.4,
'xm_h': 34.21520026085125,
'xm_v': 31.17587613521955,
'xm_w': 15.659971822785682,
'xr_h': 33.07320337042791,
'xr_v': 29.25388711043971,
'xr_w': 13.5,
'xt_h': 35.74855563619494,
'xt_v': 33.299364009371466,
'xt_w': 18.944010614572072,
'y_mlg': 2.47,
'y_n': 2.6,
'ym_h': 1.9611423076663264,
'ym_w': 5.569532204800901,
'yt_h': 4.594361473577855,
'yt_w': 14.037494434549208,
'z_lg': -2.0,
'z_n': 0.0,
'z_tailstrike': -0.84,
'zm_h': 4.42748459846653,
'zm_v': 2.070850918999471,
'zm_w': 0.4872709290626237,
'zr_h': 4.359,
'zr_v': 0.0,
'zr_w': 0.0,
'zt_h': 4.519438637980579,
'zt_v': 4.358807176281144,
'zt_w': 1.2281216273313065,
'frac_xcg_ae': None}

""" dt.balance(airplane)
print('\n##### Test Balance Module ##### \n')
print(airplane['xcg_fwd']) """

# 3.10

# dt.landing_gear(airplane)
# print('\n##### Test Landing Gear Module ##### \n')
# print(f"frac_nlg_fwd={airplane['frac_nlg_fwd']}\n")
# print(f"frac_nlg_aft={airplane['frac_nlg_aft']}\n")
# print(f"alpha_tipback={airplane['alpha_tipback']}\n")
# print(f"alpha_tailstrike={airplane['alpha_tailstrike']}\n")
# print(f"phi_overturn={airplane['phi_overturn']}\n")


# 4.3

W0_guess = 467500.00000000000000
T0_guess = 140250.00000000000000

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
    airplane['frac_xcg_ae'] = fracs
    dt.analyze(airplane, print_log=False, plot=False, W0_guess=W0_guess, T0_guess=T0_guess, review=False)  
    deltaS_wlan[i] = airplane['deltaS_wlan']
    SM_aft[i] = airplane['SM_aft']
    SM_fwd[i] = airplane['SM_fwd']
    CLv[i] = airplane['CLv']
    frac_nlg_fwd[i] = airplane['frac_nlg_fwd']
    frac_nlg_aft[i] = airplane['frac_nlg_aft']
    alpha_tipback[i] = airplane['alpha_tipback']
    alpha_tailstrike[i] = airplane['alpha_tailstrike']
    phi_overturn[i] = airplane['phi_overturn']
    T0_vecs[i] = airplane['T0vec']   

airplane['frac_xcg_ae'] = 0.423549784987987958498562
dt.analyze(airplane, print_log=True, plot=False, W0_guess=W0_guess, T0_guess=T0_guess, review=True)   

plt.figure()
plt.plot(frac_xcg_AllElses*100, deltaS_wlan)
plt.scatter(airplane['frac_xcg_ae']*100, airplane['deltaS_wlan'], marker='*', color='green', s=100)
plt.title('deltaS_wlan')
plt.axhline(0, linestyle='--', color = 'red', label='Limite')
plt.xlabel('Fração Xcg_ae [%]')
plt.legend()
    
plt.figure()
plt.plot(frac_xcg_AllElses*100, SM_aft)
plt.scatter(airplane['frac_xcg_ae']*100, airplane['SM_aft'], marker='*', color='green', s=100)
plt.title('SM_aft')
plt.axhline(0.05, linestyle='--', color = 'red', label='Limite')
plt.xlabel('Fração Xcg_ae [%]')
plt.legend()

plt.figure()
plt.plot(frac_xcg_AllElses*100, SM_fwd)
plt.title('SM_fwd')
plt.axhline(0.30, linestyle='--', color = 'red', label='Limite')
plt.scatter(airplane['frac_xcg_ae']*100, airplane['SM_fwd'], marker='*', color='green', s=100)
plt.xlabel('Fração Xcg_ae [%]')
plt.legend()

plt.figure()
plt.plot(frac_xcg_AllElses*100, CLv)
plt.title('CLv')
plt.axhline(0.75, linestyle='--', color = 'red', label='Limite')
plt.scatter(airplane['frac_xcg_ae']*100, airplane['CLv'], marker='*', color='green', s=100)
plt.xlabel('Fração Xcg_ae [%]')
plt.legend()

plt.figure()
plt.plot(frac_xcg_AllElses*100, frac_nlg_fwd)
plt.title('frac_nlg_fwd')
plt.axhline(0.18, linestyle='--', color = 'red', label='Limite')
plt.scatter(airplane['frac_xcg_ae']*100, airplane['frac_nlg_fwd'], marker='*', color='green', s=100)
plt.xlabel('Fração Xcg_ae [%]')
plt.legend()

plt.figure()
plt.plot(frac_xcg_AllElses*100, frac_nlg_aft)
plt.title('frac_nlg_aft [%]')
plt.axhline(0.05, linestyle='--', color = 'red', label='Limite')
plt.scatter(airplane['frac_xcg_ae']*100, airplane['frac_nlg_aft'], marker='*', color='green', s=100)
plt.xlabel('Fração Xcg_ae [%]')
plt.legend()


plt.figure()
plt.plot(frac_xcg_AllElses*100, np.rad2deg(alpha_tipback))
plt.scatter(airplane['frac_xcg_ae']*100, np.rad2deg(airplane['alpha_tipback']), marker='*', color='green', s=100)
plt.title('alpha_tipback')
plt.axhline(15, linestyle='--', color = 'red', label='Limite')
plt.xlabel('Fração Xcg_ae [%]')
plt.legend()


plt.figure()
plt.plot(frac_xcg_AllElses*100, np.rad2deg(alpha_tailstrike))
plt.scatter(airplane['frac_xcg_ae']*100, np.rad2deg(airplane['alpha_tailstrike']), marker='*', color='green', s=100)
plt.title('alpha_tailstrike')
plt.axhline(10, linestyle='--', color = 'red', label='Limite')
plt.xlabel('Fração Xcg_ae [%]')
plt.legend()


plt.figure()
plt.plot(frac_xcg_AllElses*100, np.rad2deg(phi_overturn))
plt.scatter(airplane['frac_xcg_ae']*100, np.rad2deg(airplane['phi_overturn']), marker='*', color='green', s=100)
plt.title('phi_overturn')
plt.axhline(63, linestyle='--', color = 'red', label='Limite')
plt.xlabel('Fração Xcg_ae [%]')
plt.legend()


plt.show()






