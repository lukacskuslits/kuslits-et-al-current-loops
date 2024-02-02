import scipy.io
import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from sklearn.preprocessing import PolynomialFeatures
from sklearn import linear_model
from sklearn.metrics import mean_squared_error as mse
import math
import numpy.matlib
import gc

# convert loop and core flow images to aitoff projection ----------------------------

import numpy as np
import matplotlib.pyplot as plt

im = plt.imread("loop_pos_2019.jpg")
im = 255-im
x = np.linspace(-np.pi,np.pi,im.shape[1])
y = np.linspace(-np.pi/2,np.pi/2,im.shape[0])
X,Y = np.meshgrid(x,y)

plt.subplot(111, projection="aitoff")
plt.pcolormesh(X,Y[::-1],im[:,:,2])
plt.show()



# single induced fields -------------------------------------------------------------

magnins_field_0 = pd.read_csv('induced_br_test_magnins.csv', skiprows=8, sep=',', header=None)
magncond_field_0 = pd.read_csv('induced_br_test_magncond.csv', skiprows=8,  sep=',', header=None)
magnins_field_0.columns = ['fi','Brad']
magncond_field_0.columns = ['fi','Brad']

plt.plot(magnins_field_0['Brad'])
plt.plot(magncond_field_0['Brad'])

# # ------------ parametric sweep moving dipole -------------------------------------
# stationary_field_dipole = pd.read_csv('mf_norm_stat.csv')
# moving_field_dipole = pd.read_csv('mf_norm_move.csv')
# plt.plot(moving_field_dipole['Bnorm'] - stationary_field_dipole['Bnorm'])

#TODO: ne futtasd, mert rosszul irja felul a npy file-okat!
# # ------------ parametric sweep stationary dipole ---------------------------------




#other file: field_param_sweep.csv
diff_field = pd.DataFrame({'fi': [], 'Brad': [], 'dtI': []})

for ii in range(15):
    print(ii)
    diff_field_part = pd.read_csv('induced_br_sweep_total'+str(ii)+'.csv', skiprows=8, sep=',', header=None)
    diff_field_part.columns = ['fi','Brad','dtI']
    diff_field = pd.concat([diff_field, diff_field_part], sort=False)
# diff_field[:5000].to_csv('induced_br_sweep_total13.csv', index=False)
# diff_field[5000:].to_csv('induced_br_sweep_total14.csv', index=False)
# diff = np.array(list(diff_field['Brad'].values[:992]) + [0,0,0,0,0,0,0,0] + list(diff_field['Brad'].values[993:]))
diff = diff_field['Brad'].values

#print(len(list(diff_field[:5000].fi.unique())))
Roc = 3.48e6
fis = np.array([0] + list(diff_field.fi.unique())) #  + [Roc*math.pi]
#print(max(list(diff_field.fi.unique())))
dtIs = np.array([10**(-3), 1])
rs = [2e5, 2.5e5, 3e5, 3.5e5, 4e5, 4.5e5, 5e5, 5.5e5, 6e5, 6.5e5, 7e5, 7.5e5, 8e5]
#rs.reverse()
rs = np.array(rs)
hs = [1e5, 3e5, 5e5, 7e5, 9e5]
#hs.reverse()
hs = np.array(hs)
params = np.meshgrid(hs, rs, dtIs, fis)
fi_vals = np.array([params[3].reshape(-1)]).T
dtI_vals = np.array([params[2].reshape(-1)]).T
h_vals = np.array([params[1].reshape(-1)]).T
r_vals = np.array([params[0].reshape(-1)]).T
print(len(diff))
print(max(diff))
params_all = np.concatenate([fi_vals, dtI_vals, h_vals, r_vals], axis=1)
print('params: ', len(params_all[:, 0]))
dIt = params_all[:, 1]

params_all = np.delete(params_all, 1, 1)
# params_all = params_all[:-7]
plt.plot(diff)
# # plt.plot(dtI_vals)
plt.show()
np.savez('diff.npz', diff)
np.savez('params_all.npz', params_all)


#---------------- parametric sweep test --------------------------------------------------------------------------------
# diff_loaded = np.load('diff_test.npz')['arr_0']
# params_loaded = np.load('params_test.npz')['arr_0']
# params_loaded = params_loaded[:2000,:]
# diff_loaded = diff_loaded[:2000]
diff_field = pd.read_csv('induced_br_test.csv', skiprows=8, sep=',', header=None)
diff_field.columns = ['fi','Brad','dtI']
#diff = np.array(list(diff_field['Brad'].values[:992]) + [0,0,0,0,0,0,0,0] + list(diff_field['Brad'].values[993:]))
diff = diff_field['Brad'].values
#print(len(list(diff_field.fi.unique())))
params_all = np.empty((0,4))
dtI_list = [10**(-3), 1, 10**(-3)]#[1, 1, 1e-3, 5e-1]
r_list = [3.25e5, 4.25e5, 5.25e5]#[1e6, 9e5, 4e5, 6e5]
h_list = [1.5e5, 3.5e5, 5.5e5]#[8e5, 8e5, 2e5, 6e5]
fis = np.array([0] + list(diff_field.fi.unique())) #  + [Roc*math.pi]
#print(max(list(diff_field.fi.unique())))
for ii in range(3):
    dtIs = np.array([dtI_list[ii]])
    rs = np.array([r_list[ii]])
    hs = np.array([h_list[ii]])
    params = np.meshgrid(hs, rs, dtIs, fis)
    fi_vals = np.array([params[3].reshape(-1)]).T
    dtI_vals = np.array([params[2].reshape(-1)]).T
    h_vals = np.array([params[1].reshape(-1)]).T
    r_vals = np.array([params[0].reshape(-1)]).T
    print(len(diff))
    print(max(diff))
    params = np.concatenate([fi_vals, dtI_vals, h_vals, r_vals], axis=1)
    params_all = np.concatenate([params_all, params], axis=0)
print('params: ', params_all[:, 0])
dIt = params_all[:, 1]
params_all = np.delete(params_all, 1, 1)
#diff = np.concatenate([diff_loaded, diff], axis=0)
diff = diff/dIt
#plt.plot(diff)
#plt.plot(params_all[:,1]/8e5)
np.savez('diff_test.npz', diff)
np.savez('params_test.npz', params_all)
#-----------------------------------------------------------------------------------------------------------------------

# plt.colorbar()
#
# diff = total_field_0['Brad']-primary_field_0['Brad']
# plt.plot(diff)
# plt.show()

#year = 3.15*10**7
# svar = 10000*10**(-9)
# dIt = svar/(30*year)
#
# # ------------ parametric sweep -----------------------
primary_field_r = pd.read_csv('primary_field_sweep_r.csv', sep=',')
primary_field_i = pd.read_csv('primary_field_sweep_i.csv', sep=',')
primary_field_h = pd.read_csv('primary_field_sweep_h.csv', sep=',')
total_field_r = pd.read_csv('total_field_sweep_r.csv', sep=',')
total_field_i = pd.read_csv('total_field_sweep_i.csv', sep=',')
total_field_h = pd.read_csv('total_field_sweep_h.csv', sep=',')

total_field = pd.concat([total_field_r, total_field_i.drop(columns=['fi','Brad']), total_field_h.drop(columns=['fi','Brad'])], axis=1)
total_field.to_csv('total_field_sweep.csv')
primary_field = pd.concat([primary_field_r, primary_field_i.drop(columns=['fi','Brad']), primary_field_h.drop(columns=['fi','Brad'])], axis=1)
primary_field.to_csv('primary_field_sweep.csv')
#------------------------------------------------------
# plt.figure(1)
# plt.plot(primary_field['Brad'])
# plt.figure(2)
# plt.plot(total_field['Brad'])
diff = total_field['Brad']-primary_field['Brad']
plt.show()

primary_field_0 = pd.read_csv('primary_field_sweep.csv', sep=',')
total_field_0 = pd.read_csv('total_field_sweep_0.csv', sep=',')

# diff = total_field_0['Brad']-primary_field_0['Brad']

plt.plot(primary_field['Brad'])
plt.plot(total_field['Brad'])
plt.plot(primary_field['r']*10**(-18), label='r')
plt.plot(primary_field['h']*10**(-18))
plt.plot(primary_field['dtI'])
plt.legend()
plt.show()


#----------------------- benchmark ----------------------------------------------------
benchmark_comsol = pd.read_csv('test_benchmark4.csv', sep=',', skiprows=8, header=None)
benchmark_comsol.columns = ['fi','Brad']
benchmark_analytic = scipy.io.loadmat('benchmark4.mat')
# plt.plot(benchmark_comsol['fi'], benchmark_analytic['res'], label='analytic solution')
# plt.plot(benchmark_comsol['fi'], benchmark_comsol['Brad'], label='numerical solution')
analytic = np.array(list(benchmark_analytic['benchmark_res']),)
analytic = analytic[:, 0]
analytic = analytic[2::1000]
numerical = np.array(list(benchmark_comsol['Brad']),)
numerical = numerical/1000
del benchmark_analytic, benchmark_comsol
gc.collect()
fig, axs = plt.subplots(1,1)
plt.rcParams["figure.figsize"] = [7.50, 3.50]
plt.rcParams["figure.autolayout"] = True
axs.plot(analytic, label='Analitic solution', linestyle='dotted', linewidth=10, color='blue')
axs.plot(numerical, label='Numerical solution', linewidth=5, color='red')
labels = [item.get_text() for item in axs.get_xticklabels()]
for i in range(len(labels)):
    labels[i] = (i-1)*36
axs.set_xticklabels(labels)
axs.tick_params(axis='x', labelsize=20)
axs.tick_params(axis='y', labelsize=20)
axs.set_title(' R=9e5 [m], I=3e12 [A], z=3e5 [m], $H_{ref}$ = $R_{CMB}+2175$ [km]', fontsize=30)
axs.set_xlabel('$\phi$ geographic distance from source [°]', fontsize=20)
axs.set_ylabel('Br($\phi$) [T]', fontsize=20)
handles1, labels1 = axs.get_legend_handles_labels()
by_label1 = dict(zip(labels1, handles1))
axs.legend(by_label1.values(), by_label1.keys(), prop={'size': 20}, loc='upper right')#, bbox_to_anchor=(0.04, 0.82, 1., .102))
plt.show()



# benchmark_comsol = pd.read_csv('benchmark_earth_like.csv', sep=',')
# benchmark_analytic = scipy.io.loadmat('Hz_anal.mat')
# plt.plot(benchmark_comsol['z'], benchmark_analytic['full_Hz'], label='analytic solution')
# plt.plot(benchmark_comsol['z'], benchmark_comsol['Hz'], label='numerical solution')
# plt.title('Comparison under Core-like conditions (sigma=5e8[S/m])')
# plt.legend()
# plt.show()
#
# benchmark_comsol = pd.read_csv('benchmark_earth_like_flip_cond.csv', sep=',')
# benchmark_analytic = scipy.io.loadmat('Hz_anal.mat')
# plt.plot(benchmark_comsol['z'], benchmark_analytic['full_Hz'], label='analytic solution')
# plt.plot(benchmark_comsol['z'], benchmark_comsol['Hz'], label='numerical solution')
# plt.title('Comparison when the mantle is supposed to be very conductive in the numerical benchmark (sigma=5e8[S/m])')
# plt.legend()
# plt.show()


#------------------- current loop vs magnetic dipole -------------------------------------------------------------------

Br_loop5 = np.array(list(scipy.io.loadmat('res_loop_1e5.mat')['res_loop']))
Br_dip5 = 4*math.pi*np.array(list(scipy.io.loadmat('res_dip_1e5.mat')['res_dip']))
Br_loop6 = np.array(list(scipy.io.loadmat('res_loop_1e6.mat')['res_loop']))
Br_dip6 = 4*math.pi*np.array(list(scipy.io.loadmat('res_dip_1e6.mat')['res_dip']))

gc.collect()
fig, axs = plt.subplots(1,2)
plt.rcParams["figure.figsize"] = [7.50, 3.50]
plt.rcParams["figure.autolayout"] = True
axs[0].plot(Br_loop5, label='Köráram R=1e5 [m] I=2.2e8 [A]', linestyle='dotted', linewidth=10, color='blue')
axs[0].plot(Br_dip5, label='Dipólus M=6.912e18 [Am$^2$]', linewidth=5, color='red')
labels = [item.get_text() for item in axs[0].get_xticklabels()]
for i in range(len(labels)):
    labels[i] = (i-1)*36
axs[0].set_xticklabels(labels)
axs[0].tick_params(axis='x', labelsize=20)
axs[0].tick_params(axis='y', labelsize=20)
axs[0].set_title('A', fontsize=30)
axs[0].set_xlabel('$\phi$ szögtávolság a forrástól [°]', fontsize=20)
axs[0].set_ylabel('Br($\phi$)', fontsize=20)
handles1, labels1 = axs[0].get_legend_handles_labels()
by_label1 = dict(zip(labels1, handles1))
axs[0].legend(by_label1.values(), by_label1.keys(), prop={'size': 20}, loc='upper right')#, bbox_to_anchor=(0.04, 0.82, 1., .102))
axs[1].plot(Br_loop6, label='Köráram R=1e6 [m] I=2.2e8 [A]', linestyle='dotted', linewidth=10, color='blue')
axs[1].plot(Br_dip6, label='Dipólus M=6.912e20 [Am$^2$]', linewidth=5, color='red')
labels = [item.get_text() for item in axs[1].get_xticklabels()]
for i in range(len(labels)):
    labels[i] = (i-1)*36
axs[1].set_xticklabels(labels)
axs[1].tick_params(axis='x', labelsize=20)
axs[1].tick_params(axis='y', labelsize=20)
axs[1].set_title('B', fontsize=30)
axs[1].set_xlabel('$\phi$ szögtávolság a forrástól [°]', fontsize=20)
axs[1].set_ylabel('Br($\phi$)', fontsize=20)
handles, labels = axs[1].get_legend_handles_labels()
by_label = dict(zip(labels, handles))
axs[1].legend(by_label.values(), by_label.keys(), prop={'size': 20}, loc='upper right')#bbox_to_anchor=(0.01, 0.82, 1., .102))
axs[1].ticklabel_format(style='sci', axis='y' ,scilimits=(0,0))
plt.show()


#--------- results of comsol simulation r3e5 z3e5 dtI5e-1---------------------------------------------------------------

cols = ['t', 'Brad']
annum=3.15e7
#induced_field_w_phi = pd.read_csv('induced_br_w_phi_time1e6.csv', sep=',', skiprows=8)
induced_field_w_time_phi_0 = pd.read_csv('induced_br_w_time_phi_0.csv', sep=',', skiprows=8)
induced_field_w_time_phi_half_pi = pd.read_csv('induced_br_w_time_phi_05pi.csv', sep=',', skiprows=8)
induced_field_w_time_phi_pi = pd.read_csv('induced_br_w_time_phi_pi.csv', sep=',', skiprows=8)
induced_field_w_time_phi_0.columns = cols
induced_field_w_time_phi_half_pi.columns = cols
induced_field_w_time_phi_pi.columns = cols
induced_field_w_time_phi_0 = induced_field_w_time_phi_0[:180]
induced_field_w_time_phi_half_pi = induced_field_w_time_phi_half_pi[:180]
induced_field_w_time_phi_pi = induced_field_w_time_phi_pi[:180]

primary_field_phi_0 = -0.015*np.ones([180,1]).squeeze()
primary_field_phi_90 = 1.38*10**(-5)*np.ones([180,1]).squeeze()
primary_field_phi_180 = 6.6*10**(-6)*np.ones([180,1]).squeeze()

gc.collect()
fig, axs = plt.subplots(1,3)
plt.rcParams["figure.figsize"] = [7.50, 3.50]
plt.rcParams["figure.autolayout"] = True
axs[0].plot(induced_field_w_time_phi_0['t']/annum, induced_field_w_time_phi_0['Brad'], linewidth=5, color='blue')
axs[0].plot(induced_field_w_time_phi_0['t']/annum, primary_field_phi_0, linewidth=5, linestyle='dotted', color='red')
axs[0].tick_params(axis='x', labelsize=20)
axs[0].tick_params(axis='y', labelsize=20)
axs[0].set_title('$\phi$=0°', fontsize=30)
axs[0].set_xlabel('t szimulációs idő [év]', fontsize=20)
axs[0].set_ylabel('Br$_{tot}$(t)-Br$_p$(t) [T]', fontsize=20)
axs[0].ticklabel_format(style='sci', axis='both', scilimits=(0, 0))
#
axs[1].plot(induced_field_w_time_phi_half_pi['t']/annum, induced_field_w_time_phi_half_pi['Brad'], linewidth=5, color='blue')
axs[1].plot(induced_field_w_time_phi_half_pi['t']/annum, primary_field_phi_90, linewidth=5, linestyle='dotted', color='red')
axs[1].tick_params(axis='x', labelsize=20)
axs[1].tick_params(axis='y', labelsize=20)
axs[1].set_title('$\phi$=90°', fontsize=30)
axs[1].set_xlabel('t szimulációs idő [év]', fontsize=20)
axs[1].ticklabel_format(style='sci', axis='both', scilimits=(0, 0))
#
axs[2].plot(induced_field_w_time_phi_pi['t']/annum, induced_field_w_time_phi_pi['Brad'], linewidth=5, color='blue')
axs[2].plot(induced_field_w_time_phi_half_pi['t']/annum, primary_field_phi_180, linewidth=5, linestyle='dotted', color='red')
axs[2].tick_params(axis='x', labelsize=20)
axs[2].tick_params(axis='y', labelsize=20)
axs[2].set_title('$\phi$=180°', fontsize=30)
axs[2].set_xlabel('t szimulációs idő [év]', fontsize=20)
axs[2].ticklabel_format(style='sci', axis='both', scilimits=(0, 0))
plt.show()

#-----------------------------------------------------------------------------------------------------------------------

#total_field_05_mu_100_t= pd.read_csv('induction_100_1_higher_mu.csv', sep=',')
primary_field_100_t = pd.read_csv('insulation_100_1.csv', sep=',')
total_field_1_mu_1e6_t = pd.read_csv('induction_1e6_1.csv', sep=',')
total_field_10_t_I_15 = pd.read_csv('induction_10_I_15.csv', sep=',').astype(float)
primary_field_10_t_I_15 = pd.read_csv('insulation_10_I_15.csv', sep=',').astype(float)
total_field_10_t_I_15_b = pd.read_csv('induction_10_I_15_b.csv', sep=',').astype(float)
primary_field_10_t_I_15_b = pd.read_csv('insulation_10_I_15_b.csv', sep=',').astype(float)

fig, axs = plt.subplots(2, 2)
# axs[0, 0].plot(total_field_1_mu_1e6_t['Bnorm']-primary_field_1e6_t['Bnorm'])
# axs[0, 0].set_title('Mu 1; t 1e6; dI/dt 1.05')
#axs[0, 1].plot(total_field_05_mu_100_t['Bnorm']-primary_field_100_t['Bnorm'])
#axs[0, 1].set_title('Mu 1.05; t 100; dI/dt 1.05')
axs[1, 0].plot(total_field_10_t_I_15['Bnorm']-primary_field_10_t_I_15['Bnorm'])
axs[1, 0].set_title('Mu=1; t=10; dI/dt=1.5; dt=0.1; sample=1000')
axs[1, 1].plot(total_field_10_t_I_15_b['Bnorm']-primary_field_10_t_I_15_b['Bnorm'])
axs[1, 1].set_title('Mu=1; t=10; dI/dt=1.5; dt=0.01; sample=10000')
plt.show()
#
# #
# total_vs_prim_1e6 = pd.read_csv('bnorm_45_25.csv', sep=',')
# total_per_prim_1e6 = pd.read_csv('bprop_45_25.csv', sep=',')
#
# fig1, axs = plt.subplots(2, 1)
# axs[0].plot(total_vs_prim_1e6['Bnorm_pri'], label='Total field')
# axs[0].plot(total_vs_prim_1e6['Bnorm_tot'], label='Primary field')
# axs[0].legend(loc='lower left')
# axs[0].set_title('Mu 1.1; t 1e6; dI/dt 1.05')
# axs[1].plot(total_per_prim_1e6['Bp/Btot'], label='Primary field/Total field')
# axs[1].legend(loc='lower left')
# axs[1].set_title('Mu 1.1; t 1e6; dI/dt 1.05')
# plt.show()


mu = 4*math.pi*10**(-7)
I = 31500*10**(-9)
R = 3*10**5

B = mu*I/(2*R)
print(B)


#------------------------------------------- synthetic results - cross correlations ------------------------------------
import numpy as np
from scipy import signal


def corr_2d(array1, array2):
    corr = np.zeros((45, 90))
    corr21 = np.zeros((45, 90))
    corr22 = np.zeros((45, 90))
    ma1 = np.mean(array1)
    ma2 = np.mean(array2)
    for mm in range(45):
        for nn in range(90):
            corr[mm,nn] = (array1[mm, nn]-ma1)*(array2[mm, nn]-ma2)
            corr21[mm,nn] = (array1[mm, nn]-ma1)**2
            corr22[mm, nn] = (array2[mm, nn] - ma2)**2
    corr_2d = sum(sum(corr))/np.sqrt(sum(sum(corr21))*sum(sum(corr22)))
    return corr_2d


indexes = [10,91,121,710]
corr_pos = np.zeros((len(indexes),1))
corr_depth = np.zeros((len(indexes),1))
corr_rad = np.zeros((len(indexes),1))
corr_I = np.zeros((len(indexes),1))
corr_dtI = np.zeros((len(indexes),1))


ii = 0
for index in indexes:
    output = np.load('output_'+str(index)+'.npy')
    truth = np.load('label_'+str(index)+'.npy')
    output_pos = np.squeeze(output[0,0,:,:])
    output_depth = np.squeeze(output[0,1,:,:])
    output_rad = np.squeeze(output[0,2,:,:])
    output_I = np.squeeze(output[0,3,:,:])
    output_dtI = np.squeeze(output[0,4,:,:])

    true_pos = np.squeeze(truth[0,:,:])
    true_depth = np.squeeze(truth[1,:,:])
    true_rad = np.squeeze(truth[2,:,:])
    true_I = np.squeeze(truth[3,:,:])
    true_dtI = np.squeeze(truth[4,:,:])

    corr_pos[ii] = corr_2d(output_pos, true_pos)
    corr_depth[ii] = corr_2d(output_depth, true_depth)
    corr_rad[ii] = corr_2d(output_rad, true_rad)
    corr_I[ii] = corr_2d(output_I, true_I)
    corr_dtI[ii] = corr_2d(output_dtI, true_dtI)
    ii = ii+1

output = np.load('output_10.npy')
truth = np.load('label_10.npy')
output_pos = np.squeeze(output[0,0,:,:])
true_pos = np.squeeze(truth[0,:,:])

plt.figure(1)
plt.imshow(output_pos)
plt.figure(2)
plt.imshow(true_pos)

output_pos = np.squeeze(output[0,4,:,:])
true_pos = np.squeeze(truth[4,:,:])

plt.figure(1)
plt.imshow(output_pos)
plt.figure(2)
plt.imshow(true_pos)



SourcesNewModel = [15,20,35,65,85]

corr_pos = np.insert(corr_pos,[1],np.mean([corr_pos[0],corr_pos[1]]))
corr_depth = np.insert(corr_depth,[1],np.mean([corr_depth[0],corr_depth[1]]))
corr_rad = np.insert(corr_rad,[1],np.mean([corr_rad[0],corr_rad[1]]))
corr_I = np.insert(corr_I,[1],np.mean([corr_I[0],corr_I[1]]))
corr_dtI = np.insert(corr_dtI,[1],np.mean([corr_dtI[0],corr_dtI[1]]))

##

import matplotlib.pyplot as plt
corr_pos = [0.8422618561066929, 0.8526829818210404, 0.8631041075353879, 0.6740997205951182, 0.6723300897553502]
corr_depth = [0.8538766382815957, 0.8610204088801974, 0.8681641794787991, 0.6903726409975961, 0.6975223971608122]
corr_rad = [0.7546442258310492, 0.7908358597667038, 0.8270274937023584, 0.6432681573159427, 0.6519590268913915]
corr_I = [0.8540662269222554, 0.859566688757154, 0.8650671505920526, 0.6793126774867799, 0.6883633995211349]
corr_dtI = [0.8220283539276287, 0.8100445643452288, 0.7980607747628288, 0.5277529259815933, 0.4751501863479285]

plt.rcParams["figure.figsize"] = [7.50, 3.50]
plt.rcParams["figure.autolayout"] = True
fig1, axs = plt.subplots(2, 1)
axs[0].plot(SourcesNewModel, corr_pos, label='Földrajzi hely $(\Phi,\Lambda)$', marker="p", markersize=20)
axs[0].plot(SourcesNewModel, corr_depth, label='Földközépponti sugártávolság $r$', marker="p", markersize=20)
axs[0].plot(SourcesNewModel, corr_rad, label='Forrás átmérő $R$', marker="p", markersize=20)
axs[0].plot(SourcesNewModel, corr_I, label='Forrás áramerősség $I$', marker="p", markersize=20)
axs[0].plot(SourcesNewModel, corr_dtI, label='Forrás áramaváltozás $dI/dt$', marker="p", markersize=20)
axs[0].tick_params(axis='x', labelsize=16)
axs[0].tick_params(axis='y', labelsize=16)
axs[0].set_ylim([0,1])
#axs[0].set_title('A', fontsize=30)
axs[0].set_xlabel('Források darabszáma a modellben', fontsize=16)
axs[0].set_ylabel('Keresztkorrelációs együttható értéke', fontsize=16)
handles1, labels1 = axs[0].get_legend_handles_labels()
by_label1 = dict(zip(labels1, handles1))
axs[0].legend(by_label1.values(), by_label1.keys(), prop={'size': 16}, loc='lower center')
plt.show()

#-------------------------cross-correlations with results from other models
import matplotlib.pyplot as plt

SourcesNewModel = [1970, 1975, 1980, 1985, 1990, 1995, 2000, 2005, 2010, 2015, 2019]
corr_pos = [-1.4, -1, 0,0,0, 0,0,0,  -6, -2, 0.8]#[0.8, -2, -6, -1, -1.4]
corr_depth = [25, 18, 0,0,0, 0,0,0,   28, 24, 26] #[26, 24, 28, 18, 25]
corr_rad = [-8, -12, 0,0,0, 0,0,0,  -14, -2, -10]#[-10, -2, -14, -12, -8]
corr_I = [13, 12, 0,0,0, 0,0,0,  15, 18, 11]#[11, 18, 15, 12, 13]

corr_field_flow_grad = [17,15.1, 0,0,0, 0,0,0, 16.6,22,14.3]


plt.rcParams["figure.figsize"] = [7.50, 3.50]
plt.rcParams["figure.autolayout"] = True
fig1, axs = plt.subplots(1, 2)
axs[0].plot(SourcesNewModel, corr_pos, label='CMB radial field magnitudes', marker="p", markersize=20)
axs[0].plot(SourcesNewModel, corr_depth, label='CMB radial field horizontal gradients', marker="p", markersize=20)
axs[0].plot(SourcesNewModel, corr_rad, label='Core surface flow magnitudes', marker="p", markersize=20)
axs[0].plot(SourcesNewModel, corr_I, label='Core surface flow horizontal gradients', marker="p", markersize=20)
axs[0].tick_params(axis='x', labelsize=16)
axs[0].tick_params(axis='y', labelsize=16)
axs[0].set_ylim([-50,50])
#axs[0].set_title('A', fontsize=30)
axs[0].set_xlabel('Year (COVOBS model epoch)', fontsize=16)
axs[0].set_ylabel('Cross correlation against inferred loop positions [%]', fontsize=16)
handles1, labels1 = axs[0].get_legend_handles_labels()
by_label1 = dict(zip(labels1, handles1))
axs[0].legend(by_label1.values(), by_label1.keys(), prop={'size': 16}, loc='lower center')
plt.show()

#------------------------------------------- results showing how many loops were discovered (Our method vs Alldredge's)-

import matplotlib.pyplot as plt

SourcesOurs = [5,10,15,20,30,50,70,90,110,130]
discovered_ours_t2 = [0.2, 0.0, 0.0, 0.05, 0.802, 0.718, 0.66, 0.628, 0.583, 0.56]
SourcesAlldredge = [5,10,15,20,30,50,70,90,110,130]
discovered_alldredge_t2 = [0.4, 0.25, 0.37, 0.15, 0.068, 0.082, 0.073, 0.059, 0.054, 0.05]

plt.rcParams["figure.figsize"] = [7.50, 3.50]
plt.rcParams["figure.autolayout"] = True
fig1, axs = plt.subplots(2, 1)
axs[0].plot(SourcesOurs, discovered_ours_t2, label='Using the DANN trained UNet', marker="p", markersize=20)
axs[0].plot(SourcesAlldredge, discovered_alldredge_t2, label='Using radial field maxima (after Alldredge (1987))', marker="p", markersize=20)
axs[0].tick_params(axis='x', labelsize=16)
axs[0].tick_params(axis='y', labelsize=16)
axs[0].set_ylim([0,1])
axs[0].set_title('Threshold value: $|\hat{\phi}-\phi|+|\hat{\Lambda}-\Lambda| < 8°$', fontsize=30)
axs[0].set_xlabel('Number of loops in the model', fontsize=16)
axs[0].set_ylabel('Ratio of "discovered" loops', fontsize=16)
handles1, labels1 = axs[0].get_legend_handles_labels()
by_label1 = dict(zip(labels1, handles1))
axs[0].legend(by_label1.values(), by_label1.keys(), prop={'size': 16}, loc='center')
plt.show()

SourcesOurs = [5,10,15,20,30,50,70,90,110,130]
discovered_ours = [0.6, 0.7, 0.6, 0.65, 0.922, 0.9, 0.909, 0.927, 0.924, 0.93]
SourcesAlldredge = [5,10,15,20,30,50,70,90,110,130]
discovered_alldredge = [1, 0.9, 0.83, 0.825, 0.607, 0.555, 0.563, 0.568, 0.569, 0.596]

plt.rcParams["figure.figsize"] = [7.50, 3.50]
plt.rcParams["figure.autolayout"] = True
axs[1].plot(SourcesOurs, discovered_ours, label='Using the DANN trained UNet', marker="p", markersize=20)
axs[1].plot(SourcesAlldredge, discovered_alldredge, label='Using radial field maxima (after Alldredge (1987))', marker="p", markersize=20)
axs[1].tick_params(axis='x', labelsize=16)
axs[1].tick_params(axis='y', labelsize=16)
axs[1].set_ylim([0,1])
axs[1].set_title('Threshold value: $|\hat{\phi}-\phi|+|\hat{\Lambda}-\Lambda| < 40°$', fontsize=30)
axs[1].set_xlabel('Number of loops in the model', fontsize=16)
axs[1].set_ylabel('Ratio of "discovered" loops', fontsize=16)
handles1, labels1 = axs[1].get_legend_handles_labels()
by_label1 = dict(zip(labels1, handles1))
axs[1].legend(by_label1.values(), by_label1.keys(), prop={'size': 16}, loc='lower center')
plt.show()


#------------------------------------------- results of final reconstruction--------------------------------------------
#SourcesFreeOrigModel = [5,15,30]
#ErrFreeOrigModel = [0.15,0.23,0.27]
SourcesRadialOrigModel = [5,15,30]
ErrRadialOrigModel = [0.14,0.15,0.16]
SourcesNewModel = [15,20,35,65,85]
ErrNewModel = [0.08,0.09,0.11,0.12,0.13]

#NRMSFreeOrigModel = [0.04,0.07,0.09]
NRMSRadialOrigModel = [0.04,0.06,0.08]
NRMSNewModel = [0.05,0.06,0.06,0.07,0.08]

import matplotlib.pyplot as plt

plt.rcParams["figure.figsize"] = [7.50, 3.50]
plt.rcParams["figure.autolayout"] = True
fig1, axs = plt.subplots(2, 1)
axs[0].plot(SourcesNewModel, ErrNewModel, label="Updated current loop models reconstructed with the refined inference algorithm", marker="p", markersize=16) #'Időben lineárisan változó radiális állású köráram-modell'
#axs[0].plot(SourcesFreeOrigModel, ErrFreeOrigModel, label="Previous (static) current loop model using arbitrarily aligned loops", marker="v", markersize=20) #'Egyenáramú tetszőleges állású köráram-modell'
axs[0].plot(SourcesRadialOrigModel, ErrRadialOrigModel, label="Static current loop models reconstructed with the previous inference algorithm", marker="o", markersize=16) #'Egyenáramú radiális állású köráram-modell'
axs[0].tick_params(axis='x', labelsize=14)
axs[0].tick_params(axis='y', labelsize=14)
axs[0].set_ylim([0,0.25])
#axs[0].set_title('A', fontsize=30)
axs[0].set_xlabel("Number of loops in the model", fontsize=14) #'Források darabszáma a modellben'
axs[0].set_ylabel("Parameter error ($RPE(\hat{Par},Par)$)", fontsize=14)
handles1, labels1 = axs[0].get_legend_handles_labels()
by_label1 = dict(zip(labels1, handles1))
axs[0].legend(by_label1.values(), by_label1.keys(), prop={'size': 14}, loc='upper center')
axs[1].plot(SourcesNewModel, NRMSNewModel, label="Updated current loop models reconstructed with the refined inference algorithm", marker="p", markersize=16) #'Időben lineárisan változó radiális állású köráram-modell'
#axs[1].plot(SourcesFreeOrigModel, NRMSFreeOrigModel, label="Previous (static) current loop model using arbitrarily aligned loops", marker="v", markersize=20) #'Egyenáramú tetszőleges állású köráram-modell'
axs[1].plot(SourcesRadialOrigModel, NRMSRadialOrigModel, label="Static current loop models reconstructed with the previous inference algorithm", marker="o", markersize=16) #'Egyenáramú radiális állású köráram-modell'
axs[1].tick_params(axis='x', labelsize=14)
axs[1].tick_params(axis='y', labelsize=14)
axs[1].set_ylim([0,0.25])
#axs[0].set_title('A', fontsize=30)
axs[1].set_xlabel("Number of loops in the model", fontsize=14) #'Források darabszáma a modellben'
axs[1].set_ylabel("Misfit measure ($NRMS(\hat{B},B)$)", fontsize=14) #'Az illeszkedés hibája ($NRMS(\hat{B},B)$)'
handles1, labels1 = axs[0].get_legend_handles_labels()
by_label1 = dict(zip(labels1, handles1))
axs[1].legend(by_label1.values(), by_label1.keys(), prop={'size': 14}, loc='upper center')
plt.show()

#------------------------------------------- results of final reconstruction using the Alldredge method --------------------------------------------
SourcesAlldredgeModel = [5,10,20]
IniGeogrErrFivePerc = [0.113,0.18,0.44]
IniGeogrErrTenPerc = [0.226,0.359,1.084]

import matplotlib.pyplot as plt

plt.rcParams["figure.figsize"] = [7.50, 3.50]
plt.rcParams["figure.autolayout"] = True
fig1, axs = plt.subplots(2, 1)
axs[0].plot(SourcesAlldredgeModel, IniGeogrErrFivePerc, label="$\delta{\\theta_0}, \delta{\phi_0} = 5\%$", marker="p", markersize=20) #'Időben lineárisan változó radiális állású köráram-modell'
axs[0].plot(SourcesAlldredgeModel, IniGeogrErrTenPerc, label="$\delta{\\theta_0}, \delta{\phi_0} = 10\%$", marker="v", markersize=20) #'Egyenáramú tetszőleges állású köráram-modell'
axs[0].tick_params(axis='x', labelsize=16)
axs[0].tick_params(axis='y', labelsize=16)
axs[0].set_ylim([0,1.25])
#axs[0].set_title('A', fontsize=30)
axs[0].set_xlabel("Number of current loops in the models reconstructed with the algorithm of Alldredge (1987)", fontsize=16) #'Források darabszáma a modellben'
axs[0].set_ylabel("Parameter error ($RPE(\hat{Par},Par)$)", fontsize=16)
handles1, labels1 = axs[0].get_legend_handles_labels()
by_label1 = dict(zip(labels1, handles1))
axs[0].legend(by_label1.values(), by_label1.keys(), prop={'size': 16}, loc='upper center')
plt.show()

#------------------------------------------- Convergence of RMSE using the Alldredge method (5 loops, delta phi; alpha = 10%) --------------------------------------------
rmse_vals = [79.269944, 79.094059, 78.933225, 78.769889, 78.606711, 78.469997, 78.337548, 78.211499, 78.098986,\
 77.998201, 77.905992, 77.851962, 77.803719, 77.772529, 77.744613, 77.719559, 77.696318, 77.674218, 77.652864, 77.631951]
iterations = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
import matplotlib.pyplot as plt

plt.rcParams["figure.figsize"] = [7.50, 3.50]
plt.rcParams["figure.autolayout"] = True
fig1, axs = plt.subplots(2, 1)
axs[0].plot(iterations, rmse_vals, marker="x", markersize=20)
axs[0].tick_params(axis='x', labelsize=16)
axs[0].tick_params(axis='y', labelsize=14)
axs[0].set_ylim([77.5,79.5])
#axs[0].set_title('A', fontsize=30)
axs[0].set_xlabel("Number of full iterations with the algorithm of Alldredge (1987)", fontsize=16)
axs[0].set_ylabel("Misfit of SHCs ($RMS$)", fontsize=16)
axs[0].set_xticks([5,10,15,20])
plt.show()



## --------------------------reconstruction w noise ----------------------------------------------------
CrossCorr25Pos = [0.9, 0.65, 0.5]
CrossCorr25Rads = [0.91, 0.65, 0.51]
CrossCorr25Is = [0.75, 0.55, 0.35]
CrossCorr25Depths = [0.92, 0.7, 0.55]
CrossCorr25dtIs = [0.71, 0.5, 0.3]

CrossCorr100Pos = [0.67, 0.25, 0.1]
CrossCorr100Rads = [0.7, 0.11, 0.03]
CrossCorr100Is = [0.65, 0.09, 0.01]
CrossCorr100Depths = [0.69, 0.05, -0.01]
CrossCorr100dtIs = [0.48, 0.1, 0.02]

SigmaVals = [0, 0.1, 0.2]


import matplotlib.pyplot as plt

plt.rcParams["figure.figsize"] = [7.50, 3.50]
plt.rcParams["figure.autolayout"] = True
fig1, axs = plt.subplots(1, 2)
axs[0].plot(SigmaVals, CrossCorr25Pos, label='Földrajzi hely $(\Phi,\Lambda)$', marker="p", markersize=16)
axs[0].plot(SigmaVals, CrossCorr25Depths, label='Földközépponti sugártávolság $r$', marker="p", markersize=16)
axs[0].plot(SigmaVals, CrossCorr25Rads, label='Forrás sugár $R$', marker="p", markersize=16)
axs[0].plot(SigmaVals, CrossCorr25Is, label='Forrás áramerősség $I$', marker="p", markersize=16)
axs[0].plot(SigmaVals, CrossCorr25dtIs, label='Forrás áramaváltozás $dI/dt$', marker="p", markersize=16)
axs[0].tick_params(axis='x', labelsize=16)
axs[0].tick_params(axis='y', labelsize=16)
axs[0].set_ylim([0,1])
axs[0].set_title('Rekonstrukció alakulása 25 forrás esetén', fontsize=16)
axs[0].set_xlabel('Zajszint ($\sigma$)', fontsize=16)
axs[0].set_ylabel('Keresztkorrelációk ($Corr(\hat{P}(\Phi,\Lambda,\sigma), P(\Phi,\Lambda))$)', fontsize=16)
handles1, labels1 = axs[0].get_legend_handles_labels()
#by_label1 = dict(zip(labels1, handles1))
#axs[0].legend(by_label1.values(), by_label1.keys(), prop={'size': 16}, loc='lower right')
##
axs[1].plot(SigmaVals, CrossCorr100Pos, label='Földrajzi hely $(\Phi,\Lambda)$', marker="p", markersize=16)
axs[1].plot(SigmaVals, CrossCorr100Depths, label='Földközépponti sugártávolság $r$', marker="p", markersize=16)
axs[1].plot(SigmaVals, CrossCorr100Rads, label='Forrás sugár $R$', marker="p", markersize=16)
axs[1].plot(SigmaVals, CrossCorr100Is, label='Forrás áramerősség $I$', marker="p", markersize=16)
axs[1].plot(SigmaVals, CrossCorr100dtIs, label='Forrás áramaváltozás $dI/dt$', marker="p", markersize=16)
axs[1].tick_params(axis='x', labelsize=16)
axs[1].tick_params(axis='y', labelsize=16)
axs[1].set_ylim([-0.05,1])
axs[1].set_title('Rekonstrukció alakulása 100 forrás esetén', fontsize=16)
axs[1].set_xlabel('Zajszint ($\sigma$)', fontsize=16)
axs[1].set_ylabel('Keresztkorrelációk ($CCC(\hat{p}(\Phi,\Lambda,\sigma), p(\Phi,\Lambda))$)', fontsize=16)
handles1, labels1 = axs[1].get_legend_handles_labels()
by_label1 = dict(zip(labels1, handles1))
axs[1].legend(by_label1.values(), by_label1.keys(), prop={'size': 16}, loc='upper right')
plt.savefig('err_graphs_noise.jpg')
plt.show()



# Maps for real data 2010's and 1970's ---------------------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt

map_2010 = np.load('./MagneticSourceDetectionData/output2010.npy')