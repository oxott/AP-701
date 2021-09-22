# Script for Skysonic Aircraft

# IMPORTS
import designTool as dt
import numpy
import pprint

# Inputs of the Geometry Module

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
 'Lb_v': 0.4475,             #! 0.475
 'Lc_h': 3.5, # calculado e ajustado com base no ERJ-145 e CRJ-200 #! was 3.5
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
 'dihedral_h': numpy.deg2rad(0),
 'dihedral_w': numpy.deg2rad(3.5),
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
 'sweep_h': numpy.deg2rad(25),
 'sweep_v': numpy.deg2rad(40.5), # alterado para encaixar melhor a EV e EH
 'sweep_w': numpy.deg2rad(21),
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

fokker_100 = {'AR_h': 4.64,
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
 'S_w': 93.5,
 'TO_flap_def': 0.3490658503988659,
 'TO_slat_def': 0.0,
 'W_crew': 4463.55,
 'W_payload': 95519.97,
 'altitude_altcruise': 4572,
 'altitude_cruise': 10668.0,
 'altitude_landing': 0.0,
 'altitude_takeoff': 0.0,
 'b_flap_b_wing': 0.6,
 'b_slat_b_wing': 0.75,
 'c_flap_c_wing': 1.2,
 'c_slat_c_wing': 1.05,
 'c_tank_c_w': 0.4,
 'clmax_w': 2.1,
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
 'xr_w': 13.5,
 'y_mlg': 2.47,
 'y_n': 2.6,
 'z_lg': -2.0,
 'z_n': 0.0,
 'z_tailstrike': -0.84,
 'zr_h': 4.359,
 'zr_v': 0.0,
 'zr_w': 0.0}
 
## SKYSONIC
# Execute the geometry function
dt.geometry(skysonic)

# Print updated dictionary
print('skysonic = ' + pprint.pformat(skysonic))

# Generate 3D plot
dt.plot3d(skysonic)

## Fokker 100
# Execute the geometry function
##dt.geometry(fokker_100)

# Print updated dictionary
#print('fokker_100 = ' + pprint.pformat(fokker_100))

# Generate 3D plot
#dt.plot3d(fokker_100)

