import designTool as dt
import numpy as np
import matplotlib.pyplot as plt
import pprint

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
'W_crew': 4463.55,
'W_payload': 95519.97,
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
'flap_type':'double slotted',
'h_ground': 10.668000000000001,
'k_exc_drag': 0.03,
'loiter_time': 2700,
'n_engines': 2,
'n_engines_under_wing': 0,
'range_altcruise': 370400.0,
'range_cruise': 2222400.0,
'rho_f': 804,
'slat_type':'slat',
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
'''
#? 3.3
altitude = 10668.0
Mach = 0.73

dt.geometry(airplane)
C = dt.engineTSFC(Mach, altitude, airplane)
print(C)


#TODO 3.4
print('3.4')
W0_guess = 467500.0 #Libras
T0_guess = 140250.0
dt.geometry(airplane)

We, xcg_e = dt.empty_weight(W0_guess, T0_guess, airplane)
# Print updated dictionary
print('airplane = ' + pprint.pformat(airplane))
print(f'We = {We}, xCG_e = {xcg_e}')
'''
#? 3.5
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
'W_allelse': 79475.0,
'W_crew': 4463.55,
'W_eng': 35075.863123799056,
'W_f': 68890.55789155893,
'W_h': 4819.756583850933,
'W_mlg': 17087.125,
'W_nlg': 3015.3749999999995,
'W_payload': 95519.97,
'W_v': 3962.4552,
'W_w': 34840.475533682125,
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
'flap_type':'double slotted',
'h_ground': 10.668000000000001,
'k_exc_drag': 0.03,
'loiter_time': 2700,
'n_engines': 2,
'n_engines_under_wing': 0,
'range_altcruise': 370400.0,
'range_cruise': 2222400.0,
'rho_f': 804,
'slat_type':'slat',
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

W0_guess = 467500.0
dt.geometry(airplane)

Wf, Mf_cruise = dt.fuel_weight(W0_guess, airplane)
print(f'Wf = {Wf}, Mf_cruise = {Mf_cruise}')

#? 3.6
dt.geometry(airplane)
W0_guess = 467500.0
T0_guess = 140250.0

W0, We, Wf, Mf_cruise, xcg_e = dt.weight(W0_guess, T0_guess, airplane)
print(f'W0 = {W0}, We = {We}, Wf = {Wf}, Mf_cruise = {Mf_cruise}, xcg_e = {xcg_e}')



#? HOMEWORK
#! 1
plt.rcParams['font.size'] = '16'
pesos = [airplane['W_payload'], airplane['W_crew'], Wf, airplane['W_w'],
         airplane['W_h'], airplane['W_v'], airplane['W_nlg'], 
         airplane['W_mlg'], airplane['W_allelse'], airplane['W_eng'],
         airplane['W_f']]
labels = ['Payload', 'Crew', 'Fuel', 'Wing', 'HT', 'VT', 'NLG', 'MLG', 'All Else', 'Engines', 'Fuselage']
colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k', 'Pink', 'purple', 'skyblue', 'lavenderblush']
fig = plt.figure()
ax = fig.gca()
ax.pie(pesos, colors = colors,autopct='%1.1f%%',pctdistance=1.1, labeldistance=1.2)
plt.legend(labels)
plt.title('Fokker 100')

#! 2
ARs  = np.arange(7, 14.01, 0.01)
W0_ARs = np.zeros(len(ARs))
W_es = np.zeros(len(ARs))
W_Fs = np.zeros(len(ARs))
for i, AR in enumerate(ARs):
    airplane['AR_w'] = AR
    dt.geometry(airplane)
    W0_ARs[i], W_es[i], W_Fs[i], Mf_cruise, xcg_e = dt.weight(W0_guess, T0_guess, airplane)

fig = plt.figure()
ax = fig.gca()
ax.plot(ARs, W0_ARs)
ax.grid()
ax.set_xlabel('Wing AR')
ax.set_ylabel('W0 [lbs]')

fig, axs = plt.subplots(2)
fig.suptitle('Empty and Fuel Weight')
axs[0].plot(ARs, W_es)
axs[0].grid()
axs[0].set_title('Empty Weight')
axs[1].plot(ARs, W_Fs)
axs[1].grid()
axs[1].set_title('Fuel Weight')
plt.show()
