# IMPORTS
import designTool as dt
import numpy as np
import pprint

# Inputs of the Geometry Module

airplane_nosso = {'AR_h': 4.19,
 'AR_v': 1.22,
 'AR_w': 6.92,
 'BPR': 3.04,
 'Cbase': None ,
 'Cht': 0.8,
 'Cvt': 0.081,
 'D_f': 2.79,
 'D_n': 1.5,
 'LD_flap_def': 0.6981317007977318 ,
 'LD_slat_def': 0.0,
 'L_f': 23.19,
 'L_n': 4.3,
 'Lb_v': 0.49,
 'Lc_h': 3.88,
 'MLW_frac': 0.9228915662650602, # maximum landing / MTOW
 'Mach_altcruise': 0.75,
 'Mach_cruise': 0.75,
 'S_w': 68.89,
 'TO_flap_def': 0.3490658503988659 ,
 'TO_slat_def': 0.0,
 'W_crew': 4463.55 ,
 'W_payload': 6820.52,
 'altitude_altcruise': 12192,
 'altitude_cruise': 12192,
 'altitude_landing': 0.0,
 'altitude_takeoff': 0.0,
 'b_flap_b_wing': 0.6,
 'b_slat_b_wing': 0.75,
 'c_flap_c_wing': 1.2,
 'c_slat_c_wing': 1.05,
 'c_tank_c_w': 0.4,
 'clmax_w': 2.1,
 'dihedral_h': 0,
 'dihedral_w': np.deg2rad(3.5),
 'distance_landing': 1800.0 ,
 'distance_takeoff': 1800.0 ,
 'eta_h': 1.0,
 'flap_type': 'double slotted',
 'h_ground': 10.668000000000001 ,
 'k_exc_drag': 0.03,
 'loiter_time': 2700,
 'n_engines': 2,
 'n_engines_under_wing': 0,
 'range_altcruise': 370400.0 ,
 'range_cruise': 2222400.0 ,
 'rho_f': 804,
 'slat_type': 'slat',
 'sweep_h': np.deg2rad(25),
 'sweep_v': np.deg2rad(37),
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
 'x_n': 23.2,
 'x_nlg': 3.6,
 'x_tailstrike': 23.68,
 'x_tank_c_w': 0.2,
 'xcg_crew': 2.5,
 'xcg_payload': 14.4,
 'xr_w': 9.54,
 'y_mlg': 2.47,
 'y_n': 2.6,
 'z_lg': -2.0,
 'z_n': 0.0,
 'z_tailstrike': -0.84,
 'zr_h': 4.359,
 'zr_v': 0.0,
 'zr_w': 0.0}


# Execute geometry function 
dt.geometry(airplane_nosso)

# Print updated dictionary
print('airplane = ' + pprint.pformat(airplane_nosso))

# Generate 3D plot
dt.plot3d(airplane_nosso)