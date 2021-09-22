# Exercicio 05

import designTool as dt
import numpy as np
import matplotlib.pyplot as plt
import pprint


# Parâmetros do Fokker 100
airplane = {'AR_h': 4.64,
            'AR_v': 1.27,
            'AR_w': 8.43,
            'BPR': 3.04,
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
            'TO_flap_def': 0.3490658503988659,
            'TO_slat_def': 0.0,
            'W_allelse': 72882.36785331425,
            'W_crew': 4463.55,
            'W_eng ': 35075.863123799056,
            'W_f': 68890.55789155893,
            'W_h': 4819.756583850933,
            'W_mlg': 15669.709088462563,
            'W_nlg': 2765.242780316923,
            'W_payload': 95519.97,
            'W_v': 3962.4552,
            'W_w': 33199.86999399066,
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
            'zt_w': 1.2281216273313065}


# Test Performance Module

W0 = 428711.39823255542433
Mf_cruise = 0.95569551000000

dt.geometry(airplane)
T0, T0vec, deltaS_wlan, CLmaxTO = dt.performance(W0, Mf_cruise, airplane)
# print('\n##### Test Performance Module #####\n')
# print(f'T0={T0}\n')
# print(f'T0vec={T0vec}\n')
# print(f'deltaS_wlan={deltaS_wlan}\n')
# print(f'S_wlan={S_w_lan}')
# print(f'CLmaxTO={CLmaxTO}\n')


# Test Thrust Matching Module

W0_guess = 467500.00000000000000
T0_guess = 140250.00000000000000

dt.thrust_matching(W0_guess, T0_guess, airplane)
# print('\n##### Test Thrust Matching Module ##### \n')
# print('airplane = ' + pprint.pformat(airplane))


# HomeWork
print('\n##### Homework - Fokker 100 ##### \n')
S_w_range = np.arange(80, 140.1, 0.1)
T0vec_range = np.zeros((8, len(S_w_range)))

S_w_original = airplane['S_w']

for i, S_w in enumerate(S_w_range):
    airplane['S_w'] = S_w
    dt.geometry(airplane)
    dt.thrust_matching(W0_guess, T0_guess, airplane)
    T0vec_range[:, i] = airplane['T0vec']
# print(T0vec_range)


T0_range = ['T0_TO', 'T0_cruise', 'T0_FAR_25_111', 'T0_FAR_25_121a',
            'T0_FAR_25_121b', 'T0_FAR_25_121c', 'T0_FAR_25_119', 'T0_FAR_25_121d']
fig = plt.figure()
ax = fig.gca()
for i, T in enumerate(T0_range):
    plt.plot(S_w_range, T0vec_range[i, :], label=f'{T}')
plt.axvline(S_w_original - deltaS_wlan, label='T0_LD', color='k')
plt.legend()
plt.scatter(93.5, 123215, label='Fokker 100')
plt.xlim(60, 140)
plt.xlabel('S_w [m²]')
plt.ylabel('T0 [N]')
plt.title('Diagrama de Projeto Fokker-100')
# plt.show()
# print(f'Fokker:{S_w_original - deltaS_wlan}')

# PARTE 2

# Parametros aeronave Skysonic
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
            'Lb_v': 0.475,
            'Lc_h': 3.5,  # calculado e ajustado com base no ERJ-145 e CRJ-200
            'MLW_frac': 0.9228915662650602,
            'Mach_altcruise': 0.4,  # não sabemos o que é
            'Mach_cruise': 0.75,
            'S_w': 68.89,
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
            'c_tank_c_w': 0.4,
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
            'x_mlg': 17.8,
            'x_n': 16.76,
            'x_nlg': 3.6,
            'x_tailstrike': 23.68,
            'x_tank_c_w': 0.2,
            'xcg_crew': 2.5,
            'xcg_payload': 14.4,
            'xr_w': 9.54,  # estimado com base em dados do CRJ-200
            'y_mlg': 2.47,
            'y_n': 2.6,
            'z_lg': -2.0,
            'z_n': 0.0,
            'z_tailstrike': -0.84,
            'zr_h': 3.784,
            'zr_v': 0.0,
            'zr_w': 0.0}

dt.geometry(skysonic)

print('\n##### Homework - Skysonic ##### \n')

W0 = 428711.39823255542433
Mf_cruise = 0.95569551000000

T0, T0vec, deltaS_wlan_skysonic, CLmaxTO = dt.performance(
    W0, Mf_cruise, skysonic)
W0_guess_Skysonic = 467500.00000000000000
T0_guess_Skysonic = 140250.00000000000000

S_w_range = np.arange(60, 100.1, 0.1)
T0vec_range = np.zeros((8, len(S_w_range)))

S_w_original_skysonic = skysonic['S_w']

for i, S_w in enumerate(S_w_range):
    skysonic['S_w'] = S_w
    dt.geometry(skysonic)
    dt.thrust_matching(W0_guess_Skysonic, T0_guess_Skysonic, skysonic)
    T0vec_range[:, i] = skysonic['T0vec']
# print(T0vec_range)


T0_range = ['T0_TO', 'T0_cruise', 'T0_FAR_25_111', 'T0_FAR_25_121a',
            'T0_FAR_25_121b', 'T0_FAR_25_121c', 'T0_FAR_25_119', 'T0_FAR_25_121d']
fig = plt.figure()
ax = fig.gca()
for i, T in enumerate(T0_range):
    plt.plot(S_w_range, T0vec_range[i, :], label=f'{T}')
plt.axvline(S_w_original_skysonic - deltaS_wlan_skysonic,
            label='T0_LD', color='k')
plt.scatter(51.2, 77977.3, label='ERJ-145')
plt.scatter(79.1, 82025.18, label='CRJ-200')
plt.scatter(76.4, 87629.93, label='Fokker F28')
plt.scatter(73, 106000, label='Skysonic (PW812d)', color='k', marker='*')
plt.legend()
plt.xlabel('S_w [m²]')
plt.ylabel('T0 [N]')
plt.title('Diagrama de Projeto Skysonic')
plt.show()
#print(f'Skisonic:{S_w_original_skysonic - deltaS_wlan_skysonic}')
