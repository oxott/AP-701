import designTool as dt
import numpy as np
import matplotlib.pyplot as plt
import pprint

plt.rcParams['font.size'] = '16'

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

Mach = 0.30000000000000
altitude = 10.668
n_engines_failed = 1.00
flap_def = 0.34906585039887
slat_def = 0.0
lg_down = 1
h_ground = 10.668
W0_guess = 467500.0
CD0, K, CLmax = dt.aerodynamics(Mach, altitude, n_engines_failed, flap_def, slat_def,
                    lg_down, h_ground, W0_guess, airplane, method=2)
#! AVIÃO TOP

altitude = 11000
n_engines_failed = 0
flap_def = 0.0
slat_def = 0.0
lg_down = 0
h_ground = 0.0
W0_guess = 41500*9.81

airplane['sweep_w'] = np.deg2rad(10)

Machs  = np.arange(0.6, 0.801, 0.001)
Wing_Sweep = np.arange(10, 45, 5)
CD0    = np.zeros((len(Wing_Sweep), len(Machs)))
K      = np.zeros((len(Wing_Sweep), len(Machs)))
CLmax  = np.zeros((len(Wing_Sweep), len(Machs)))


for i, sweep in enumerate(Wing_Sweep):
    airplane['sweep_w'] = np.deg2rad(sweep)
    dt.geometry(airplane)
    for j,mach in enumerate(Machs):
        CD0[i,j], K[i, j], CLmax[i, j] = dt.aerodynamics(mach, altitude, n_engines_failed, flap_def, slat_def,
                    lg_down, h_ground, W0_guess, airplane, method=2)

#! Questão 1
fig = plt.figure()
ax = fig.gca()
for i, sweep in enumerate(Wing_Sweep):
    ax.plot(Machs, CD0[i,:], label=f'{sweep}°')
ax.legend()
ax.set_xlabel('Mach')
ax.set_ylabel('CD0')
ax.set_xlim(0.6, 0.8)
ax.set_ylim(np.amin(CD0), np.amax(CD0))
ax.grid()
plt.title('Fokker 100')
#plt.show()

#! Questão 2
fig2 = plt.figure()
ax2 = fig2.gca()
ax2.plot(Wing_Sweep, CLmax[:,1])
ax2.set_xlabel('Wing Sweep [°]')
ax2.set_ylabel('CLmax')
ax2.set_xlim(10, 40)
ax2.set_ylim(np.amin(CLmax), np.amax(CLmax))
ax2.grid()
plt.legend()
plt.title('Fokker 100')
#plt.show()

#! Questão 4
airplane['sweep_w'] = 0.3045599544730105
dt.geometry(airplane)

#? Cruise
Mach = 0.75
altitude = 11000
n_engines_failed = 0
flap_def = 0.0
slat_def = 0.0
lg_down = 0
h_ground = 0
W0_guess = 41500*9.81

CD0_Cruise, K_Cruise, CLmax_Cruise = dt.aerodynamics(Mach, altitude, n_engines_failed, flap_def, slat_def,
                 lg_down, h_ground, W0_guess, airplane, method=2)

#? Takeoff
Mach = 0.2
altitude = 0
n_engines_failed = 0
flap_def = 20*np.pi/180
slat_def = 0.0
lg_down = 1
h_ground = 10.67
W0_guess = 41500*9.81

CD0_TO, K_TO, CLmax_TO = dt.aerodynamics(Mach, altitude, n_engines_failed, flap_def, slat_def,
                 lg_down, h_ground, W0_guess, airplane, method=2)

#? Landing
Mach = 0.2
altitude = 0
n_engines_failed = 0
flap_def = 40*np.pi/180
slat_def = 0.0
lg_down = 1
h_ground = 10.67
W0_guess = 38300*9.81

CD0_landing, K_landing, CLmax_landing = dt.aerodynamics(Mach, altitude, n_engines_failed, flap_def, slat_def,
                 lg_down, h_ground, W0_guess, airplane, method=2)


CLs_cruise = np.arange(-0.5, CLmax_Cruise, 0.001)
CLs_TO     = np.arange(-0.5, CLmax_TO, 0.001)
CLs_landing= np.arange(-0.5, CLmax_landing, 0.001)
CD_cruise = CD0_Cruise + K_Cruise*CLs_cruise**2
CD_TO = CD0_TO + K_TO*CLs_TO**2
CD_landing = CD0_landing + K_landing*CLs_landing**2

fig3 = plt.figure()
ax3 = fig3.gca()
ax3.plot(CD_cruise, CLs_cruise, label='Cruise')
ax3.plot(CD_landing, CLs_landing, label='Landing')
ax3.plot(CD_TO, CLs_TO, label='Takeoff')
ax3.set_xlabel('CD')
ax3.set_ylabel('CL')
ax3.grid()
ax3.legend()
ax3.set_xlim(0,0.5)
ax3.set_ylim(-0.5, 3)
plt.title('Fokker 100')
#plt.show()

#! Questão 5

#? Takeoff
Mach = 0.2
altitude = 0
n_engines_failed = 0
flap_def = 20*np.pi/180
slat_def = 0.0
lg_down = 1
h_ground = 10.67
W0_guess = 41500*9.81

CD0_5 = np.zeros(5)
K_5 = np.zeros(5)
CLmax_5 = np.zeros(5)

#? Plain Flap
airplane['flap_type'] = 'plain'
flap_def = 20*np.pi/180
airplane['max_flap_def'] = 60*np.pi/180
dt.geometry(airplane)
CD0_5[0], K_5[0], CLmax_5[0] = dt.aerodynamics(Mach, altitude, n_engines_failed, flap_def, slat_def,
                 lg_down, h_ground, W0_guess, airplane, method=2)

#? Single slotted flap:
airplane['flap_type'] = 'slotted'
flap_def = 20*np.pi/180
airplane['max_flap_def'] = 40*np.pi/180
dt.geometry(airplane)
CD0_5[1], K_5[1], CLmax_5[1] = dt.aerodynamics(Mach, altitude, n_engines_failed, flap_def, slat_def,
                 lg_down, h_ground, W0_guess, airplane, method=2)

#? Fowler flap:
airplane['flap_type'] = 'fowler'
flap_def = 15*np.pi/180
airplane['max_flap_def'] = 40*np.pi/180
dt.geometry(airplane)
CD0_5[2], K_5[2], CLmax_5[2] = dt.aerodynamics(Mach, altitude, n_engines_failed, flap_def, slat_def,
                 lg_down, h_ground, W0_guess, airplane, method=2)

#? Double slotted flap:
airplane['flap_type'] = 'double slotted'
flap_def = 20*np.pi/180
airplane['max_flap_def'] = 50*np.pi/180
dt.geometry(airplane)
CD0_5[3], K_5[3], CLmax_5[3] = dt.aerodynamics(Mach, altitude, n_engines_failed, flap_def, slat_def,
                 lg_down, h_ground, W0_guess, airplane, method=2)

#? Triple slotted flap:
airplane['flap_type'] = 'triple slotted'
flap_def = 20*np.pi/180
airplane['max_flap_def'] = 40*np.pi/180
dt.geometry(airplane)
CD0_5[4], K_5[4], CLmax_5[4] = dt.aerodynamics(Mach, altitude, n_engines_failed, flap_def, slat_def,
                 lg_down, h_ground, W0_guess, airplane, method=2)

flaps = ['Plain', 'Single Slotted', 'Fowler', 'Double Slotted','Triple Slotted']

threshold = 2.1
fig, ax = plt.subplots()
barlist = ax.bar(flaps, CLmax_5)
barlist[0].set_color('r')
barlist[1].set_color('g')
barlist[2].set_color('g')
barlist[3].set_color('g')
barlist[4].set_color('g')

ax.axhline(y=threshold,linewidth=1, color='k', label='CL = 2.1')
plt.title('Fokker 100')
plt.legend()
#plt.show()

plt.figure()
plt.title('Fokker 100')
plt.bar(flaps, CD0_5)
#plt.show()

#!NOSSO AVIÃO

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
 'Lc_h': 3.5, # calculado e ajustado com base no ERJ-145 e CRJ-200
 'MLW_frac': 0.9228915662650602,
 'Mach_altcruise': 0.4, # não sabemos o que é
 'Mach_cruise': 0.75,
 'S_w': 68.89,
 'TO_flap_def': 0.3490658503988659,
 'TO_slat_def': 0.0,
 'W_crew': 4463.55,
 'W_payload': 6820.52,
 'altitude_altcruise': 4572, #não sabemos o que é
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
 'cm_h': 2.107457619636192 ,
'cm_v': 3.4576757510555542 ,
'cm_w': 3.756317488774531 ,
'cr_h': 2.849393124273043 ,
'cr_v': 3.944978890651773 ,
'cr_w': 5.3933059334262 ,
'ct_h': 1.1112633184664868 ,
'ct_v': 2.919284379082312 ,
'ct_w': 1.267426894355157 ,
#!
 'dihedral_h': np.deg2rad(0),
 'dihedral_w': np.deg2rad(3.5),
 'distance_landing': 1800.0,
 'distance_takeoff': 1800.0,
 'eta_h': 1.0,
 'flap_type': 'double slotted',
 'h_ground': 10.668000000000001, # distância do solo para que seja computado o efeito de solo
 'k_exc_drag': 0.03,
 'loiter_time': 45*60 + ((dt.nm2m * 200) / (0.75 * 343)),
 'n_engines': 2,
 'n_engines_under_wing': 0,
 'range_altcruise': 370400.0,
 'range_cruise': dt.nm2m * 2000,
 'rho_f': 804,
 'slat_type': 'slat',
 'sweep_h': np.deg2rad(25),
 'sweep_v': np.deg2rad(40.5), # alterado para encaixar melhor a EV e EH
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
 'xr_w': 9.54, # estimado com base em dados do CRJ-200
 'y_mlg': 2.47,
 'y_n': 2.6,
 'z_lg': -2.0,
 'z_n': 0.0,
 'z_tailstrike': -0.84,
 'zr_h': 3.784,
 'zr_v': 0.0,
 'zr_w': 0.0}
dt.geometry(skysonic)

altitude = 11000
n_engines_failed = 0
flap_def = 0.0
slat_def = 0.0
lg_down = 0
h_ground = 0.0
W0_guess = 41500*9.81

Machs  = np.arange(0.6, 0.801, 0.001)
Wing_Sweep = np.arange(10, 45, 5)
CD0    = np.zeros((len(Wing_Sweep), len(Machs)))
K      = np.zeros((len(Wing_Sweep), len(Machs)))
CLmax  = np.zeros((len(Wing_Sweep), len(Machs)))

for i, sweep in enumerate(Wing_Sweep):
    skysonic['sweep_w'] = np.deg2rad(sweep)
    dt.geometry(skysonic)
    for j,mach in enumerate(Machs):
        CD0[i,j], K[i, j], CLmax[i, j] = dt.aerodynamics(mach, altitude, n_engines_failed, flap_def, slat_def,
                    lg_down, h_ground, W0_guess, skysonic, method=2)

#! Questão 1
fig = plt.figure()
ax = fig.gca()
for i, sweep in enumerate(Wing_Sweep):
    ax.plot(Machs, CD0[i,:], label=f'{sweep}°')
ax.legend()
ax.set_xlabel('Mach')
ax.set_ylabel('CD0')
ax.set_xlim(0.6, 0.8)
ax.set_ylim(np.amin(CD0), np.amax(CD0))
ax.grid()
plt.title('Skysonic')
#plt.show()

#! Questão 2
fig2 = plt.figure()
ax2 = fig2.gca()
ax2.plot(Wing_Sweep, CLmax[:,1])
ax2.set_xlabel('Wing Sweep [°]')
ax2.set_ylabel('CLmax')
ax2.grid()
ax2.set_xlim(10, 40)
ax2.set_ylim(np.amin(CLmax), np.amax(CLmax))
plt.title('Skysonic')
#plt.show()

#! Questão 4
skysonic['sweep_w'] = np.deg2rad(21)
dt.geometry(skysonic)

#? Cruise
Mach = 0.75
altitude = 11000
n_engines_failed = 0
flap_def = 0.0
slat_def = 0.0
lg_down = 0
h_ground = 0
W0_guess = 41500*9.81

CD0_Cruise, K_Cruise, CLmax_Cruise = dt.aerodynamics(Mach, altitude, n_engines_failed, flap_def, slat_def,
                 lg_down, h_ground, W0_guess, skysonic, method=2)

#? Takeoff
Mach = 0.2
altitude = 0
n_engines_failed = 0
flap_def = 20*np.pi/180
slat_def = 0.0
lg_down = 1
h_ground = 10.67
W0_guess = 41500*9.81

CD0_TO, K_TO, CLmax_TO = dt.aerodynamics(Mach, altitude, n_engines_failed, flap_def, slat_def,
                 lg_down, h_ground, W0_guess, skysonic, method=2)

#? Landing
Mach = 0.2
altitude = 0
n_engines_failed = 0
flap_def = 40*np.pi/180
slat_def = 0.0
lg_down = 1
h_ground = 10.67
W0_guess = 38300*9.81

CD0_landing, K_landing, CLmax_landing = dt.aerodynamics(Mach, altitude, n_engines_failed, flap_def, slat_def,
                 lg_down, h_ground, W0_guess, skysonic, method=2)


CLs_cruise = np.arange(-0.5, CLmax_Cruise, 0.001)
CLs_TO     = np.arange(-0.5, CLmax_TO, 0.001)
CLs_landing= np.arange(-0.5, CLmax_landing, 0.001)
CD_cruise = CD0_Cruise + K_Cruise*CLs_cruise**2
CD_TO = CD0_TO + K_TO*CLs_TO**2
CD_landing = CD0_landing + K_landing*CLs_landing**2

fig3 = plt.figure()
ax3 = fig3.gca()
ax3.plot(CD_cruise, CLs_cruise, label='Cruise')
ax3.plot(CD_landing, CLs_landing, label='Landing')
ax3.plot(CD_TO, CLs_TO, label='Takeoff')
ax3.set_xlabel('CD')
ax3.set_ylabel('CL')
ax3.grid()
ax3.legend()
ax3.set_xlim(0,0.5)
ax3.set_ylim(-0.5, 3)
plt.title('Skysonic')
#plt.show()

#! Questão 5

#? Takeoff
Mach = 0.2
altitude = 0
n_engines_failed = 0
flap_def = 20*np.pi/180
slat_def = 0.0
lg_down = 1
h_ground = 10.67
W0_guess = 41500*9.81

CD0_5 = np.zeros(5)
K_5 = np.zeros(5)
CLmax_5 = np.zeros(5)

#? Plain Flap
skysonic['flap_type'] = 'plain'
flap_def = 20*np.pi/180
skysonic['max_flap_def'] = 60*np.pi/180
dt.geometry(skysonic)
CD0_5[0], K_5[0], CLmax_5[0] = dt.aerodynamics(Mach, altitude, n_engines_failed, flap_def, slat_def,
                 lg_down, h_ground, W0_guess, skysonic, method=2)

#? Single slotted flap:
skysonic['flap_type'] = 'slotted'
flap_def = 20*np.pi/180
skysonic['max_flap_def'] = 40*np.pi/180
dt.geometry(skysonic)
CD0_5[1], K_5[1], CLmax_5[1] = dt.aerodynamics(Mach, altitude, n_engines_failed, flap_def, slat_def,
                 lg_down, h_ground, W0_guess, skysonic, method=2)

#? Fowler flap:
skysonic['flap_type'] = 'fowler'
flap_def = 15*np.pi/180
skysonic['max_flap_def'] = 40*np.pi/180
dt.geometry(skysonic)
CD0_5[2], K_5[2], CLmax_5[2] = dt.aerodynamics(Mach, altitude, n_engines_failed, flap_def, slat_def,
                 lg_down, h_ground, W0_guess, skysonic, method=2)

#? Double slotted flap:
skysonic['flap_type'] = 'double slotted'
flap_def = 20*np.pi/180
skysonic['max_flap_def'] = 50*np.pi/180
dt.geometry(skysonic)
CD0_5[3], K_5[3], CLmax_5[3] = dt.aerodynamics(Mach, altitude, n_engines_failed, flap_def, slat_def,
                 lg_down, h_ground, W0_guess, skysonic, method=2)

#? Triple slotted flap:
skysonic['flap_type'] = 'triple slotted'
flap_def = 20*np.pi/180
skysonic['max_flap_def'] = 40*np.pi/180
dt.geometry(skysonic)
CD0_5[4], K_5[4], CLmax_5[4] = dt.aerodynamics(Mach, altitude, n_engines_failed, flap_def, slat_def,
                 lg_down, h_ground, W0_guess, skysonic, method=2)

flaps = ['Plain', 'Single Slotted', 'Fowler', 'Double Slotted','Triple Slotted']

threshold = 2.1
above_threshold = np.maximum(CLmax_5 - threshold, 0)
below_threshold = np.minimum(CLmax_5, threshold)

fig, ax = plt.subplots()
barlist = ax.bar(flaps, CLmax_5)
barlist[0].set_color('r')
barlist[1].set_color('g')
barlist[2].set_color('g')
barlist[3].set_color('g')
barlist[4].set_color('g')

ax.axhline(y=threshold,linewidth=1, color='k', label='CL = 2.1')
plt.legend()
plt.title('Skysonic')
#plt.show()

plt.figure()
plt.bar(flaps, CD0_5)
plt.title('Skysonic')
plt.show()