import designTool as dt
import numpy as np
import matplotlib.pyplot as plt
import pprint

airplane = {'AR_h': 4.64 ,
'AR_v': 1.27 ,
'AR_w': 8.43 ,
'BPR': 3.04 ,
'Cbase': None ,
'Cht': 0.94 ,
'Cvt': 0.088 ,
'D_f': 3.3 ,
'D_n': 1.5 ,
'LD_flap_def': 0.6981317007977318 ,
'LD_slat_def': 0.0 ,
'L_f': 32.5 ,
'L_n': 4.3 ,
'Lb_v': 0.55 ,
'Lc_h': 4.83 ,
'MLW_frac': 0.9228915662650602 ,
'Mach_altcruise': 0.4 ,
'Mach_cruise': 0.73 ,
'S_h': 18.196687370600415 ,
'S_v': 14.959999999999999 ,
'S_w': 93.5 ,
'TO_flap_def': 0.3490658503988659 ,
'TO_slat_def': 0.0 ,
'W_crew': 4463.55 ,
'W_payload': 95519.97 ,
'altitude_altcruise': 4572 ,
'altitude_cruise': 10668.0 ,
'altitude_landing': 0.0 ,
'altitude_takeoff': 0.0 ,
'b_flap_b_wing': 0.6 ,
'b_h': 9.18872294715571 ,
'b_slat_b_wing': 0.75 ,
'b_v': 4.358807176281144 ,
'b_w': 28.074988869098416 ,
'c_flap_c_wing': 1.2 ,
'c_slat_c_wing': 1.05 ,
'c_tank_c_w': 0.4 ,
'clmax_w': 2.1 ,
'cm_h': 2.107457619636192 ,
'cm_v': 3.4576757510555542 ,
'cm_w': 3.756317488774531 ,
'cr_h': 2.849393124273043 ,
'cr_v': 3.944978890651773 ,
'cr_w': 5.3933059334262 ,
'ct_h': 1.1112633184664868 ,
'ct_v': 2.919284379082312 ,
'ct_w': 1.267426894355157 ,
'dihedral_h': 0.03490658503988659 ,
'dihedral_w': 0.08726646259971647 ,
'distance_landing': 1800.0 ,
'distance_takeoff': 1800.0 ,
'eta_h': 1.0 ,
'flap_type':'double slotted',
'h_ground': 10.668000000000001 ,
'k_exc_drag': 0.03 ,
'loiter_time': 2700 ,
'n_engines': 2,
'n_engines_under_wing': 0,
'range_altcruise': 370400.0 ,
'range_cruise': 2222400.0 ,
'rho_f': 804 ,
'slat_type':'slat',
'sweep_h': 0.4537856055185257 ,
'sweep_v': 0.715584993317675 ,
'sweep_w': 0.3045599544730105 ,
'taper_h': 0.39 ,
'taper_v': 0.74 ,
'taper_w': 0.235 ,
'tcr_h': 0.1 ,
'tcr_v': 0.1 ,
'tcr_w': 0.123 ,
'tct_h': 0.1 ,
'tct_v': 0.1 ,
'tct_w': 0.096 ,
'x_mlg': 17.8 ,
'x_n': 23.2 ,
'x_nlg': 3.6 ,
'x_tailstrike': 23.68 ,
'x_tank_c_w': 0.2 ,
'xcg_crew': 2.5 ,
'xcg_payload': 14.4 ,
'xm_h': 34.21520026085125 ,
'xm_v': 31.17587613521955 ,
'xm_w': 15.659971822785682 ,
'xr_h': 33.07320337042791 ,
'xr_v': 29.25388711043971 ,
'xr_w': 13.5 ,
'xt_h': 35.74855563619494 ,
'xt_v': 33.299364009371466 ,
'xt_w': 18.944010614572072 ,
'y_mlg': 2.47 ,
'y_n': 2.6 ,
'ym_h': 1.9611423076663264 ,
'ym_w': 5.569532204800901 ,
'yt_h': 4.594361473577855 ,
'yt_w': 14.037494434549208 ,
'z_lg': -2.0,
'z_n': 0.0 ,
'z_tailstrike': -0.84 ,
'zm_h': 4.42748459846653 ,
'zm_v': 2.070850918999471 ,
'zm_w': 0.4872709290626237 ,
'zr_h': 4.359 ,
'zr_v': 0.0 ,
'zr_w': 0.0 ,
'zt_h': 4.519438637980579 ,
'zt_v': 4.358807176281144 ,
'zt_w': 1.2281216273313065}

dt.geometry(airplane)
#print('airplane = ' + pprint.pformat(airplane))
#dt.plot3d(airplane)

Mach = 0.30000000000000
altitude = 11000
n_engines_failed = 0
flap_def = 0.0
slat_def = 0.0
lg_down = 0
h_ground = 0.0
W0_guess = 41500*9.81
airplane['sweep_w'] = np.deg2rad(10)

Machs  = np.arange(0.6, 0.8, 0.01)
Wing_Sweep = np.arange(10, 40, 5)
CD0    = np.zeros((len(Wing_Sweep), len(Machs)))
K      = np.zeros((len(Wing_Sweep), len(Machs)))
CLmax  = np.zeros((len(Wing_Sweep), len(Machs)))

for i, sweep in enumerate(Wing_Sweep):
    airplane['sweep_w'] = np.deg2rad(sweep)
    dt.geometry(airplane)
    for j,mach in enumerate(Machs):
        CD0[i,j], K[i, j], CLmax[i, j] = dt.aerodynamics(mach, altitude, n_engines_failed, flap_def, slat_def,
                    lg_down, h_ground, W0_guess, airplane, method=2)


fig = plt.figure()
ax = fig.gca()
for i, sweep in enumerate(Wing_Sweep):
    ax.plot(Machs, CD0[i,:], label=f'{sweep}°')
ax.legend()
ax.set_xlabel('Mach')
ax.set_ylabel('CD0')


fig2 = plt.figure()
ax2 = fig2.gca()
for i, sweep in enumerate(Wing_Sweep):
    ax2.plot(Machs, CLmax[i,:], label=f'{sweep}°')
ax2.legend()
ax2.set_xlabel('Mach')
ax2.set_ylabel('CLmax')
plt.show()