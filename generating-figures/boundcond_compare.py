import pandas as pd

magnind_field_tang1 = pd.read_csv('induced_br_test_boudcond_tang.csv', skiprows=8, sep=',', header=None)
magnind_field_rad1 = pd.read_csv('induced_br_test_boudcond_rad.csv', skiprows=8,  sep=',', header=None)
magnind_field_tang1.columns = ['fi','Brad']
magnind_field_rad1.columns = ['fi','Brad']
magnind_field_rad1 = magnind_field_rad1.iloc[::2,:]
magnind_field_rad1 = magnind_field_rad1.reset_index()

magnind_field_tang2 = pd.read_csv('induced_br_test_boudcond_tang_3e5.csv', skiprows=8, sep=',', header=None)
magnind_field_rad2 = pd.read_csv('induced_br_test_boudcond_rad_3e5.csv', skiprows=8,  sep=',', header=None)
magnind_field_tang2.columns = ['fi','Brad']
magnind_field_rad2.columns = ['fi','Brad']
magnind_field_rad2 = magnind_field_rad2.iloc[::2,:]
magnind_field_rad2 = magnind_field_rad2.reset_index()

magnind_field_tang3 = pd.read_csv('induced_br_test_boudcond_tang_9e5.csv', skiprows=8, sep=',', header=None)
magnind_field_rad3 = pd.read_csv('induced_br_test_boudcond_rad_9e5.csv', skiprows=8,  sep=',', header=None)
magnind_field_tang3.columns = ['fi','Brad']
magnind_field_rad3.columns = ['fi','Brad']
magnind_field_rad3 = magnind_field_rad3.iloc[::2,:]
magnind_field_rad3 = magnind_field_rad3.reset_index()

very_small1 = max(abs(magnind_field_rad1.Brad))*1e-2
very_small2 = max(abs(magnind_field_rad2.Brad))*1e-2
very_small3 = max(abs(magnind_field_rad3.Brad))*1e-2

#import math
import numpy as np
bound_effect_comp1 = np.abs(1-(very_small1+abs(magnind_field_tang1.Brad))/(very_small1 + abs(magnind_field_rad1.Brad)))
bound_effect_comp2 = np.abs(1-(very_small2+abs(magnind_field_tang2.Brad))/(very_small2 + abs(magnind_field_rad2.Brad)))
bound_effect_comp3 = np.abs(1-(very_small3+abs(magnind_field_tang3.Brad))/(very_small3 + abs(magnind_field_rad3.Brad)))
degs = np.arange(0,360,360/500)

import matplotlib.pyplot as plt
plt.plot(degs, bound_effect_comp2, linewidth=3, label='$R=3e5$')
plt.plot(degs, bound_effect_comp1, linewidth=3, label='$R=7e5$')
plt.plot(degs, bound_effect_comp3, linewidth=3, label='$R=9e5$')
plt.xlabel('Distane from the loop center at the CMB [$\circ$]', fontsize=30)
plt.ylabel('$|1-Br(\phi\')_{perfins}/Br(\phi\')_{perfcond}|$', fontsize=30)
plt.xticks(fontsize=26)
plt.yticks(fontsize=20)
plt.legend(loc='upper right', fontsize=26)
plt.show()
