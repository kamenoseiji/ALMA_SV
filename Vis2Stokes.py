execfile(SCR_DIR + 'interferometry.py')
from scipy.constants import constants
from scipy.interpolate import UnivariateSpline
from scipy.interpolate import griddata
antList = np.load(prefix + '.Ant.npy')
AzEl = np.load(prefix + '.Azel.npy')
Vis  = np.load(prefix + '.Vis.npy')

#-------- Least-Square fit for polarizatino parameters (Q, U, XYphase, Dx, Dy)
timeNum = Vis.shape[1]
PA2 = 2.0*AzEl[3]
P = np.zeros([5, 2* timeNum])
solution = np.array([0.1, 0.1, np.angle( np.mean(Vis[1])), 0.0, 0.0])   # Initial parameters
#-------- Iteration loop
for index in range(10):
    modelVis = np.r_[ np.cos(solution[2])* (-np.sin(PA2)* solution[0] + np.cos(PA2)* solution[1]) + solution[3], np.sin(solution[2])* (-np.sin(PA2)* solution[0] + np.cos(PA2)* solution[1]) + solution[4] ]
    #-------- Partial matrix
    P[0] = np.r_[-np.sin(PA2)* np.cos(solution[2]), -np.sin(PA2)* np.sin(solution[2])]
    P[1] = np.r_[ np.cos(PA2)* np.cos(solution[2]),  np.cos(PA2)* np.sin(solution[2])]
    P[2] = np.r_[ -np.sin(solution[2])* (-np.sin(PA2)* solution[0] + np.cos(PA2)* solution[1]), np.cos(solution[2])* (-np.sin(PA2)* solution[0] + np.cos(PA2)* solution[1])]
    P[3] = np.r_[np.ones([timeNum]), np.zeros([timeNum])]
    P[4] = np.r_[np.zeros([timeNum]), np.ones([timeNum])]
    #PTP = np.dot(P, P.T)
    #PTP_inv = scipy.linalg.inv(PTP)
    PTP_inv = scipy.linalg.inv(np.dot(P, P.T))
    vecVis = np.r_[ Vis[1].real, Vis[1].imag ]
    #residual = visvec - np.r_[ np.cos(solution[2])* (-np.sin(PA2)* solution[0] + np.cos(PA2)* solution[1]) + solution[3], np.sin(solution[2])* (-np.sin(PA2)* solution[0] + np.cos(PA2)* solution[1]) + solution[4] ]
    residual = vecVis - modelVis
    correction = np.dot( PTP_inv, np.dot (P, residual))
    solution   = solution + correction
#

plt.plot(AzEl[3], Vis[1].real, 'b.')
plt.plot(AzEl[3], Vis[1].imag, 'g.')
PArange = np.arange(min(AzEl[3]), max(AzEl[3]), 0.01)
plt.plot(PArange, np.cos(solution[2])* (-np.sin(2.0*PArange)* solution[0] + np.cos(2.0* PArange)* solution[1]) + solution[3], 'b-')
plt.plot(PArange, np.sin(solution[2])* (-np.sin(2.0*PArange)* solution[0] + np.cos(2.0* PArange)* solution[1]) + solution[4], 'g-')
plt.xlabel('PA [rad]'); plt.ylabel('XY (real and imaginary)'); plt.title(prefix)
text_sd = 'Q/I = %6.3f   U/I = %6.3f   XY_phase = %6.3f rad (RefAnt : %s)' % (solution[0], solution[1], solution[2], antList[0]); plt.text(min(AzEl[3]), min(Vis[1].real), text_sd, size='x-small')
np.save( prefix + 'QUXY.npy', solution )
