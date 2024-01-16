import math
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import PolynomialFeatures, SplineTransformer
from sklearn.linear_model import Ridge, Lasso
from sklearn import linear_model
from sklearn.metrics import mean_squared_error as mse
from sklearn.pipeline import make_pipeline


params_all = np.load('params_all.npz')['arr_0']
diff = np.load('diff.npz')['arr_0']
params_all = params_all[:,:]
diff = diff[:]
print(len(diff))
scale = 1/params_all.max(axis=0)
# ----------- param fit - polynomial -----------------------------------------------
# https://stackoverflow.com/questions/3433486/how-to-do-exponential-and-logarithmic-curve-fitting-in-python-i-found-only-poly
poly = PolynomialFeatures(degree=11)
param_selection = params_all
diff_selection = diff#[1000:]
w = 5e-3 + np.sqrt(abs(diff_selection))
diff_selection = np.log(2.5*10**(-3)-diff_selection)
# plt.plot(diff)
# plt.show()
# param_selection = np.concatenate([param_selection], axis=0)
# diff_selection = np.concatenate([diff_selection], axis=0)
param_selection = param_selection*scale
X_fit = poly.fit_transform(param_selection)

vector = np.array([diff_selection]).T
#generate the regression object
clf = linear_model.Ridge(alpha=10**(-7))
#preform the actual regression
clf.fit(X_fit, vector, sample_weight=w)

true = diff
preds = clf.predict(X_fit)
preds = 2.5*10**(-3) - np.exp(preds)
print(mse(preds, true)/(max(true)-min(true)))
plt.figure(1)
plt.plot(preds, color='red')
plt.scatter(list(range(len(true))), true)
plt.plot(param_selection[:, 0])
plt.plot(param_selection[:, 1])
plt.plot(param_selection[:, 2])
plt.show()
del true, preds

# ------------------------- param fit - test -------------------------------------------------------------------------

param_test = np.load('params_test.npz')['arr_0']
diff_test = np.load('diff_test.npz')['arr_0']
param_test = param_test[:,:]
diff_test = diff_test[:]
prediction_loc = param_test*scale
prediction_loc_polynomial = poly.fit_transform(prediction_loc)
preds = clf.predict(prediction_loc_polynomial)
print(mse(2.5*10**(-3) - np.exp(preds), diff_test)/(max(diff_test)-min(diff_test)))
import matplotlib
fig, ax = plt.subplots()
fig.canvas.draw()
plt.plot(2.5*10**(-3)-np.exp(preds), color='red', label='Approximate value') #'Közelítő érték'
plt.scatter(list(range(len(diff_test))), diff_test, label='Simulation result') #'Szimulációból számított érték'
plt.xlim([-1, 1500])
start, end = ax.get_xlim()
#labels = [item.get_text() for item in ax.get_xticklabels()]
labels = [0, 30, 60, 120, 0, 30, 60, 120, 0, 30, 60, 120]
# multiple = 0
# for i in range(len(labels)):
#     if (i+1) % 3 == 0:
#         multiple += 1m
#     labels[i] = i*90-180*multiple
ax.set_xticks(np.arange(start, end, 125)[:-1])
#plt.locator_params(axis='x', nbins=12)
ax.set_xticklabels(labels)
# plt.plot(prediction_loc[:,0])
# plt.plot(prediction_loc[:,1])
# plt.plot(prediction_loc[:,2])
plt.xlabel('$\phi$\' angular distance from the source [°]') #'$\phi$\' szögtávolság a forrástól [°]'
plt.ylabel('C($\phi$\') = Bi($\phi$\')*(dI/dt)^(-1)')

#plt.xticks(labels)
plt.legend()
matplotlib.rcParams.update({'font.size': 32})
plt.show()
import scipy.io
scipy.io.savemat('coef.mat', {'coef': clf.coef_})
scipy.io.savemat('powers.mat', {'powers': poly.powers_})
scipy.io.savemat('weights.mat', {'weights': w})


# # ------------------------- param fit - spline -------------------------------------------------------------------------
# params_all = np.load('params_all.npz')['arr_0']
# diff = np.load('diff.npz')['arr_0']
# scale = 1/params_all.max(axis=0)
# # ----------- param fit - polynomial -----------------------------------------------
# # https://stackoverflow.com/questions/3433486/how-to-do-exponential-and-logarithmic-curve-fitting-in-python-i-found-only-poly
# param_selection = params_all
# w = 5e-3 + np.sqrt(abs(diff))
# diff = np.log(1.5*10**(-3)-diff)
# pipeline = make_pipeline(SplineTransformer(n_knots=12, degree=11), Ridge(alpha=1e-7))
# pipeline.fit(param_selection, diff, ridge__sample_weight=w)
# preds = pipeline.predict(param_selection)
# plt.plot(preds, color='red')
# plt.scatter(list(range(len(diff))), diff)
# plt.show()

# preds = 1.5*10**(-3) - np.exp(preds)
# diff = 1.5*10**(-3) - np.exp(diff)
# plt.scatter(list(range(len(diff))), diff)
# plt.plot(preds, color='red')
# plt.show()
