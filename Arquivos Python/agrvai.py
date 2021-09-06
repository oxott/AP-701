import numpy as np
kg2lbs = 2.20462
ft2m = 0.3048
m2ft = (1/ft2m)
m2toft2 = m2ft**2
Nz = 1.5*2.5                                        #? Load Factor
S_w = 93.5
S_csw = 0.15*(S_w*m2toft2)
W0_guess = 467500.0             #Libras
AR_w = 8.43
tcr_w = 0.123
taper_w = 0.235
sweep_w = 0.3045599544730105

W_w = 0.0051*pow(W0_guess*Nz, 0.557)*pow(S_w*m2toft2, 0.649)*\
    pow(AR_w, 0.55)*pow(tcr_w, -0.4)*pow(1 + taper_w, 0.1)*\
        pow(np.cos(sweep_w), -1)*pow(S_csw, 0.1)
        
print(W_w)