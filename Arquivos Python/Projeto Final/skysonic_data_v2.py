import numpy as np
import designTool as dt
# DADOS DA AERONAVE SKYSONIC

# Atualização: 23/09/2021
# Atualizado por: Felipe Miotto


skysonic = {'S_w': 68.89,  # Wing area [m2]
            'AR_w': 7.7,  # Wing aspect ratio
            'taper_w': 0.35,  # Wing taper ratio
            'sweep_w': np.deg2rad(28),  # Wing sweep [rad]
            'dihedral_w': np.deg2rad(3.5),  # Wing dihedral [rad]
            # Longitudinal position of the wing (with respect to the fuselage nose) [m]
            'xr_w': 9.54,
            # Vertical position of the wing (with respect to the fuselage nose) [m]
            'zr_w': 0.0,
            'tcr_w': 0.12,  # t/c of the root section of the wing
            'tct_w': 0.12,  # t/c of the tip section of the wing

            'clmax_w': 2.3134,  # Maximum lift coefficient of wing airfoil

            'TO_flap_def': 20*np.pi/180,  # Takeoff flap deflection [rad]
            'LD_flap_def': 40*np.pi/180,  # Landing flap deflection [rad]
            'flap_type': 'triple slotted',  # Flap type
            'c_flap_c_wing': 1.2,  # chord_with_deflected_flaps/chord_with_retracted_flaps
            'b_flap_b_wing': 0.6,  # Fraction of the wing span occupied by flaps

            'TO_slat_def': 0*np.pi/180,  # Takeoff slat deflection [rad]
            'LD_slat_def': 0*np.pi/180,  # Landing slat deflection [rad]
            'slat_type': 'slat',  # Slat type
            'c_slat_c_wing': 1.05,  # chord_with_deflected_slats/chord_with_retracted_slats
            'b_slat_b_wing': 0.75,  # Fraction of the wing span occupied by slats


            'Cht': 0.8,  # Horizontal tail volume coefficient
            # Non-dimensional lever of the horizontal tail (lever/wing_mac)
            'Lc_h': 3.7,
            'AR_h': 4.19,  # HT aspect ratio
            'taper_h': 0.52,  # HT taper ratio
            'sweep_h': np.deg2rad(30),  # HT sweep [rad]
            'dihedral_h': np.deg2rad(0),  # HT dihedral [rad]
            'zr_h': 4.289,  # Vertical position of the HT [m]
            'tcr_h': 0.1,  # t/c of the root section of the HT
            'tct_h': 0.1,  # t/c of the tip section of the HT
            'eta_h': 1.0,  # Dynamic pressure factor of the HT

            'Cvt': 0.081,  # Vertical tail volume coefficient
            # Non-dimensional lever of the vertical tail (lever/wing_span)
            'Lb_v': 0.38,
            'AR_v': 1.22,  # VT aspect ratio
            'taper_v': 0.7,  # VT taper ratio
            'sweep_v': np.deg2rad(43),  # VT sweep [rad]
            'zr_v': 0.0,  # Vertical position of the VT [m]
            'tcr_v': 0.1,  # t/c of the root section of the VT
            'tct_v': 0.1,  # t/c of the tip section of the VT

            'L_f': 23.19,  # Fuselage length [m]
            'D_f': 2.79,  # Fuselage diameter [m]

            # Longitudinal position of the nacelle frontal face [m]
            'x_n': 16.76,
            'y_n': 2.6,  # Lateral position of the nacelle centerline [m]
            'z_n': 0.0,  # Vertical position of the nacelle centerline [m]
            'L_n': 4.3,  # Nacelle length [m]
            'D_n': 1.5,  # Nacelle diameter [m]

            'n_engines': 2,  # Number of engines
            'n_engines_under_wing': 0,  # Number of engines installed under the wing
            'BPR': 5.5,  # Engine bypass ratio
            # Base engine TSFC [1/s] (use 'None' for Howe's values)
            'Cbase': 0.5,

            # Longitudinal position of the nose landing gear [m]
            'x_nlg': 3.6,
            # Longitudinal position of the main landing gear [m]
            'x_mlg': 17.8,
            'y_mlg': 2.47,  # Lateral position of the main landing gear [m]
            'z_lg': -2.0,  # Vertical position of the landing gear [m]
            # Longitudinal position of critical tailstrike point [m]
            'x_tailstrike': 23.68,
            # Vertical position of critical tailstrike point [m]
            'z_tailstrike': -0.84,

            'c_tank_c_w': 0.515,  # Fraction of the wing chord occupied by the fuel tank
            'x_tank_c_w': 0.2,  # Fraction of the wing chord where fuel tank starts


            # Distance to the ground for ground effect computation [m]
            'h_ground': 35.0*dt.ft2m,
            'k_exc_drag': 0.03,  # Excrescence drag factor

            # Altitude for takeoff computation [m]
            'altitude_takeoff': 0.0,
            'distance_takeoff': 1800.0,  # Required takeoff distance [m]

            # Altitude for landing computation [m]
            'altitude_landing': 0.0,
            # Required landing distance [m] (The actual Fokker100 distance is 1350 m but it is very restrictive compared to the historical regression. Therefore I kept the same TO distance since the aircraft should takeoff and land at the same runway)
            'distance_landing': 1800.0,
            'MLW_frac': 38300/41500,  # Max Landing Weight / Max Takeoff Weight

            'altitude_cruise': 12192,  # Cruise altitude [m]
            'Mach_cruise': 0.75,  # Cruise Mach number
            'range_cruise': dt.nm2m * 2000,  # Cruise range [m]

            # Loiter time [s]
            'loiter_time': 45*60 + ((dt.nm2m * 200) / (0.75 * 343)),

            'altitude_altcruise': 4572,  # Alternative cruise altitude [m]
            'Mach_altcruise': 0.4,  # Alternative cruise Mach number
            'range_altcruise': 370400.0,  # Alternative cruise range [m]

            'W_payload': 15036.67*dt.lb2N,  # Payload weight [N]
            # Longitudinal position of the Payload center of gravity [m]
            'xcg_payload': 14.4,

            'W_crew': 1240*dt.lb2N,  # Crew weight [N]
            # Longitudinal position of the Crew center of gravity [m]
            'xcg_crew': 2.5,

            'rho_f': 804,  # Fuel density kg/m3 (This is Jet A-1)
            }


# Parâmetros que precisam ser definidos para o skysonic:

# 'D_n': 1.5 -> Diametro da Nacelle
# 'L_n': 4.3,  # Nacelle length [m]
# 'x_nlg': 3.6, Longitudinal position of the nose landing gear [m]
# 'x_mlg': 17.8, # Longitudinal position of the main landing gear [m]
# 'y_mlg': 2.47,  # Lateral position of the main landing gear [m]
# 'z_lg': -2.0,  # Vertical position of the landing gear [m]
# 'x_tailstrike': 23.68, # Longitudinal position of critical tailstrike point [m]
# 'z_tailstrike': -0.84, # Vertical position of critical tailstrike point [m]
# 'c_tank_c_w': 0.515,  # Fraction of the wing chord occupied by the fuel tank
# 'x_tank_c_w': 0.2,  # Fraction of the wing chord where fuel tank starts
# 'xcg_crew': 2.5, # Longitudinal position of the Crew center of gravity [m]
# 'xcg_payload': 14.4, # Longitudinal position of the Payload center of gravity [m]

# Obs: Verificar o angulo de enflechamento do EH. Se deve ser maior ou menor que o da asa

# Obs: Sempre rodar os módulos de geometry e aerodinâmica antes de outras análises
