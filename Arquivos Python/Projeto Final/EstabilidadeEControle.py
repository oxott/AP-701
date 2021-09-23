import designTool as dt
import numpy as np
import matplotlib.pyplot as plt
from skysonic_data_v2 import skysonic

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