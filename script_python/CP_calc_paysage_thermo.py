import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from matplotlib import cm, ticker
import time
# alpha = 0.1406250000e-3
# beta  = 115.2
# L     = 0.05
# h     = 4 * np.sqrt( alpha / (2 * beta)) # epaisseur interface
# Nx = int(5 * L /h) ; Ny = Nx
# print('epaisseur interface :' , float(h))
# print('nbre maille necessaire :' , float(Nx))



# x0 = 0.0
# y0 = 0.012
# R  = 0.09486832980505137
# xlim_l = -0.025 ; xlim_r = -xlim_l
# ylim_b = -0.025 ; ylim_t = -ylim_b
# def z_func(x,y,str_type):
	# if str_type == 'bulle':
		# return (-0.5*np.tanh((np.sqrt((x-x0)**2+(y+y0)**2)-R**2)/h))
	# elif str_type == 'base':
		# return (0.5*np.tanh(640.*(y+0.0125-0.01666666667*np.exp(-0.2304e6*(x+0.025)*(x+0.025)))))
	# elif str_type == 'vague':
		# return (-0.5*np.tanh((y+0.0150-0.005*np.cos(40*np.pi*x))/(h*np.sqrt(2))))
	# else :
		# print('erreur')
		

# x = np.linspace(xlim_l,xlim_r,Nx)
# y = np.linspace(ylim_b,ylim_t,Ny)
# X,Y = np.meshgrid(x, y) # grid of point
# Z = z_func(X, Y,'vague') # evaluation of the function on the grid
# extent = np.min(x), np.max(x), np.min(y), np.max(y)
# im = plt.imshow(Z,extent = extent,cmap='jet') # drawing the function
# plt.hlines(y, xmin = xlim_l, xmax = xlim_r , linewidth=0.5 ) 
# plt.vlines(x, ymin = ylim_b, ymax = ylim_t , linewidth=0.5 ) 

# plt.colorbar(im) # adding the colobar on the right
# plt.title('')

# plt.show()
# plt.close()

plot_pot_chim = 'non'
use_eigen_val = 'non'
use_det_hess = 'oui'


def deg2rad(_Value_deg):
	return (_Value_deg * np.pi / 180)
	
def rotation(theta):
	Mat_Rot = sp.Matrix([ [ sp.cos(theta), sp.sin(theta)],
						  [-sp.sin(theta),  sp.cos(theta)] ])
	return Mat_Rot

def construct_paraboloid(_coord, a0, a1):
	return ( _coord.row(0) / a0)**2 + ( _coord.row(1) / a1)**2

theta =  0
x_u = sp.Symbol("x_u")
x_zr = sp.Symbol("x_zr")

x_u_eq_ox = 0. #eq 1
x_zr_eq_ox = 1  #eq 1
x_u_eq_met =  0.06822576740446558	  #0.0233103  #eq 0
x_zr_eq_met = 0#0 #eq 0
a_oxU  = 1       #demi puit phase oxyde U
a_oxZr = 1  #demi puit phase oxyde Zr
a_metU  = 0.5#0.5    #demi puit phase met U
a_metZr = 0.8#2.8   4   #demi puit phase met Zr

# InitialPoints
xUInit0 = 0.4121405750798722 #0.41#0.68#0.15#0.68#0.15
xUInit1 = 0
xZrInit0 = 1 - xUInit0 #0.59#0.32#0.85#0.32#0.85
xZrInit1 = 0

# Equilibrium points
xUEq0= x_u_eq_ox      #0 #0.5 xMiscible a l'equilibre 0 = goutte
xZrEq0 = x_zr_eq_ox   #1 # 0.1 xImmiscible a l'equilibre 0 = goutte
xUEq1 = x_u_eq_met # 0.038085# 0.051 #0.00939712169# 0.2 #0.68 xMiscible a l'equilibre 1 = phase continue
xZrEq1 = x_zr_eq_met # 0. # 0.7 xImmiscible a l'equilibre 1 = phase continue

X_ox = sp.Matrix([[ x_u - x_u_eq_ox, x_zr - x_zr_eq_ox]])
X_ox = sp.transpose(X_ox)
res_ox = rotation(deg2rad(0)) * X_ox 
X_met = sp.Matrix([[ x_u - x_u_eq_met, x_zr - x_zr_eq_met]])
X_met = sp.transpose(X_met)
res_met = rotation(deg2rad(0)) * X_met 




P_ox  = ((1 / a_oxU ) * res_ox.row(0))**2 + ((1 / a_oxZr)  * res_ox.row(1))**2			#PARABOLOIDE
P_met = ((1 / a_metU) * res_met.row(0))**2 + ((1 / a_metZr) * res_met.row(1))**2		#PARABOLOIDE

x = np.linspace(0,1,100)
result = np.zeros((x.size, x.size), dtype=np.float64)									#Omega -> potentiel analytique

G = sp.Function('f')(x_u,x_zr)
G = P_met[0,0] * P_ox[0,0]
t1 = time.process_time()
for i, x1 in enumerate(x):
	for j, x2 in enumerate(x):
		if x1+x2<1:
			result[i][j] = G.subs({x_u:x[i], x_zr:x[j]}).evalf()
		if result[i][j] < 1e-10:
			result[i][j] = np.NaN
		
t2 = time.process_time()
print('temps necessaire pour tracer le paysage : ',t2-t1,'s')

X, Y = np.meshgrid(x, x)
Z = np.ma.masked_invalid(result)
Z = np.transpose(result)  ## transpose necessary for python management of lists of lists

G = sp.Function('f')(x_u,x_zr)
G = P_met[0,0] * P_ox[0,0]


t1 = time.process_time()
Hess_Gibbs = sp.matrices.dense.hessian(G, (x_u,x_zr))
if use_eigen_val == 'oui':
# Calcul les valeurs propres de la matrice hessienne --> zone spinodale comprise dans la zone instable
	eigenVal_formel = Hess_Gibbs.eigenvals(error_when_incomplete = True, multiple = True)
	xUeigenValNeg =[]
	xZreigenValNeg = []
	for i, x1 in enumerate(x):
		for j, x2 in enumerate(x):
			if x1+x2<1:
				eigenVal1 = eigenVal_formel[0].subs({x_u:x[i], x_zr:x[j]}).evalf()
				eigenVal2 = eigenVal_formel[1].subs({x_u:x[i], x_zr:x[j]}).evalf()
				if eigenVal1 < 0 or eigenVal2 < 0:
					print("xU = ",x[i], "xZr = ",x[j]," eigen",eigenVal1, eigenVal2)
					xUeigenValNeg=np.append(xUeigenValNeg,x[i])
					xZreigenValNeg=np.append(xZreigenValNeg,x[j]) 


if use_det_hess == 'oui':
# determine la zone instable telle que det(d2G/dx1dx2) < 0 (det(mat_Hessienne) < 0) 
	abs_hess_formel = (Hess_Gibbs)
	det_abs_hess_formel = abs_hess_formel.det()
	xUValNeg = [] ; xZrValNeg =[]
	for i, x1 in enumerate(x):
		for j, x2 in enumerate(x):
			if x1+x2<1:
				val_det_hess = det_abs_hess_formel.subs({x_u:x[i], x_zr:x[j]}).evalf()
				if val_det_hess < 0:
					# print("xU = ",x[i], "xZr = ",x[j]," eigen",val_det_hess)
					xUValNeg=np.append(xUValNeg,x[i])
					xZrValNeg=np.append(xZrValNeg,x[j])

t2 = time.process_time()


print('temps de calcul pour la determination de la zone instable : ',t2-t1, 's')
size = 9
plt.rcParams['figure.figsize'] = (size * (np.max(X) - np.min(X)) / (np.max(Y) - np.min(Y)), size)
fig, axes = plt.subplots(1, 1, constrained_layout=True)
ax = axes
SMALL_SIZE = 8
MEDIUM_SIZE = 10
BIGGER_SIZE = 20  # 50
plt.rcParams.update({'font.size': BIGGER_SIZE})
plt.rcParams['text.usetex'] = True
plt.rc('axes', titlesize=BIGGER_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize=BIGGER_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=BIGGER_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize=BIGGER_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize=BIGGER_SIZE)  # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
ax.grid(True)
locator = ticker.LogLocator(subs=range(1, 10))
# locator = ticker.LinearLocator(numticks=100)
plotd = plt.contourf(X, Y, Z, cmap='jet', locator=locator)
plt.rcParams['text.usetex'] = True
if use_eigen_val == 'oui': plt.plot(xUeigenValNeg, xZreigenValNeg, marker = "+",color = 'red',linestyle = 'None')
if use_det_hess == 'oui': plt.plot(xUValNeg, xZrValNeg,marker =  "x",color ='blue',linestyle = 'None')
plt.plot([x_u_eq_met, x_u_eq_ox], [x_zr_eq_met, x_zr_eq_ox], 'b-', lw=2)  # blue straight line
plt.plot([xUInit0, xUInit1], [xZrInit0, xZrInit1], 'r-', lw=2)  # Red straight line
plt.ylabel(r"$\phi_{Zr}$",fontsize=BIGGER_SIZE)
plt.xlabel(r"$\phi_U$",fontsize=BIGGER_SIZE)
plt.xticks(fontsize=BIGGER_SIZE)
plt.yticks(fontsize=BIGGER_SIZE)
cb_ax = fig.add_axes([0.8, 0.3, 0.02, 0.6])
cbar = fig.colorbar(plotd, format='%4.3e', cax=cb_ax)
plt.savefig('landscape_spinodal.eps')
plt.legend
plt.show()
plt.close()











#####################
# POTENTIEL CHIMIQUE 
#####################


if plot_pot_chim == 'oui':
	poten_chim_ox  = (P_ox * sp.diff(P_met,x_u ) + P_met * sp.diff(P_ox,x_u ))[0,0]
	poten_chim_met = (P_ox * sp.diff(P_met,x_zr) + P_met * sp.diff(P_ox,x_zr))[0,0]
	for i, x1 in enumerate(x):
		for j, x2 in enumerate(x):
			if x1+x2<1:
				result[i][j] = poten_chim_ox.subs({x_u:x[i], x_zr:x[j]}).evalf()
			if result[i][j] < 1e-10:
				result[i][j] = np.NaN
	plt.rcParams['figure.figsize'] = (size * (np.max(X) - np.min(X)) / (np.max(Y) - np.min(Y)), size)
	fig, axes = plt.subplots(1, 1, constrained_layout=True)
	ax = axes
	plt.rcParams.update({'font.size': BIGGER_SIZE})
	plt.rcParams['text.usetex'] = True
	plt.rc('axes', titlesize=BIGGER_SIZE)  # fontsize of the axes title
	plt.rc('axes', labelsize=BIGGER_SIZE)  # fontsize of the x and y labels
	plt.rc('xtick', labelsize=BIGGER_SIZE)  # fontsize of the tick labels
	plt.rc('ytick', labelsize=BIGGER_SIZE)  # fontsize of the tick labels
	plt.rc('legend', fontsize=BIGGER_SIZE)  # legend fontsize
	plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
	ax.grid(True)
	locator = ticker.LogLocator(subs=range(1, 10))
	plotd = plt.contourf(X, Y, result, cmap='jet', locator=locator)
	plt.rcParams['text.usetex'] = True
	plt.ylabel(r"$\phi_{Zr}$",fontsize=BIGGER_SIZE)
	plt.xlabel(r"$\phi_U$",fontsize=BIGGER_SIZE)
	plt.xticks(fontsize=BIGGER_SIZE)
	plt.yticks(fontsize=BIGGER_SIZE)
	plt.title(r"$\mu_{U}$")
	cb_ax = fig.add_axes([0.8, 0.3, 0.02, 0.6])
	cbar = fig.colorbar(plotd, format='%4.3e', cax=cb_ax)
	plt.savefig('mu_ox.eps')
	plt.legend
	plt.show()
	plt.close()



	for i, x1 in enumerate(x):
		for j, x2 in enumerate(x):
			if x1+x2<1:
				result[i][j] = poten_chim_met.subs({x_u:x[i], x_zr:x[j]}).evalf()
			if result[i][j] < 1e-10:
				result[i][j] = np.NaN
	plt.rcParams['figure.figsize'] = (size * (np.max(X) - np.min(X)) / (np.max(Y) - np.min(Y)), size)
	fig, axes = plt.subplots(1, 1, constrained_layout=True)
	ax = axes
	plt.rcParams.update({'font.size': BIGGER_SIZE})
	plt.rcParams['text.usetex'] = True
	plt.rc('axes', titlesize=BIGGER_SIZE)  # fontsize of the axes title
	plt.rc('axes', labelsize=BIGGER_SIZE)  # fontsize of the x and y labels
	plt.rc('xtick', labelsize=BIGGER_SIZE)  # fontsize of the tick labels
	plt.rc('ytick', labelsize=BIGGER_SIZE)  # fontsize of the tick labels
	plt.rc('legend', fontsize=BIGGER_SIZE)  # legend fontsize
	plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
	ax.grid(True)
	locator = ticker.LogLocator(subs=range(1, 10))
	plotd = plt.contourf(X, Y, result, cmap='jet', locator=locator)
	plt.rcParams['text.usetex'] = True
	plt.ylabel(r"$\phi_{Zr}$",fontsize=BIGGER_SIZE)
	plt.xlabel(r"$\phi_U$",fontsize=BIGGER_SIZE)
	plt.xticks(fontsize=BIGGER_SIZE)
	plt.yticks(fontsize=BIGGER_SIZE)
	plt.title(r"$\mu_{Zr}$")
	cb_ax = fig.add_axes([0.8, 0.3, 0.02, 0.6])
	cbar = fig.colorbar(plotd, format='%4.3e', cax=cb_ax)
	plt.savefig('mu_ox.eps')
	plt.legend
	plt.show()
	plt.close()
