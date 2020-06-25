# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 13:26:48 2019

@author: B. Bonnet
"""

#faudra ensuite : - trouver phi critique, passer à phi compression au dessus de phi critique


import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy import misc
import time
from mpl_toolkits.mplot3d import Axes3D

start = time.time()


#FUNCTIONS
#------------------------------------------------------------------------------


def mu(phi): # (eq.18)
    return mu_b*( 1 + (1.25*phi)/ (1- phi/phi_cp) )**2


def dPI_dphi(phi):
    return misc.derivative(PI,phi,dx=1e-6)
def dPI0_dphi(phi):
    return misc.derivative(PI0,phi,dx=1e-6)

def H(phi):  # Happel function (eq.6)
    return (6 + 4*phi**(5/3.) ) / ( 6 - 9*phi**(1./3.) + 9*phi**(5/3.) - 6*phi**2 )


def D(phi): # (eq.6)
    return Vp*dPI_dphi(phi) / (6*np.pi*mu(phi)*a*H(phi))
def D_0(phi): # (eq.6)
    return Vp*dPI0_dphi(phi) / (6*np.pi*mu(phi)*a*H(phi))


def G_s(phi_w): #simplified expression of G(phib,phiw) (D and mu: constants) (s stands for simplified) (eq.21)
    return (D0**2/mu_b) * ( np.log(phi_w)*(phi_w - phi_b) - phi_w*(np.log(phi_w)-1) + phi_b*(np.log(phi_b)-1) )
    # return ((D(phi_w)**2)/mu(phi_w))*( np.log(phi_w)*(phi_w - phi_b) - phi_w*(np.log(phi_w)-1) + phi_b*(np.log(phi_b)-1) )
    
    
def f1(phi): # function of inner integral (eq.21)
    return D(phi) / (phi*mu(phi))


def int1(phi0, phi_wi): # inner integral (eq.21)
    MAX= 50
    int0= 0
    phi= np.linspace(phi0,phi_wi,MAX)
    for i in range(MAX-1):
        int0 = int0 +  ( f1(phi[i]) + f1(phi[i+1]) ) * (phi[i+1]-phi[i]) / 2. 
    return int0


def f2(phi, phi_wi):
    return D(phi)*int1(phi, phi_wi) # function of outter integral (eq.21)


def G(phi_wi): # (eq21) (double integral solved by trapezoidal integration)
    MAX= 20 # 20 seems a good balance between processing time and and precision
    int2= 0
    phi= np.linspace(phi_b,phi_wi,MAX)
    for i in range(MAX-1):
        int2= int2 + ( f2(phi[i], phi_wi)+f2(phi[i+1], phi_wi) ) * (phi[i+1]-phi[i]) / 2.
    return int2


def eq20(Q): #(eq.20)
    K= rho*f*GG/( phi_b*np.pi*(R**3)*(Vwi**2) )  # GG = G(phi_wi)
    return Q0 - Q - (Q**2)*K # develloped expression of the line below
    #return Vwi*np.sqrt(Q0-Q) - (Q/R)*np.sqrt( rho*f*G_s(phi_wi) / (phi_b*np.pi*R) )
    
def inttrapz(i): # Only the last integration is performed, the other ones are kept in memory in intz (list)
    int_temp= (f3(phi_w[i]) + f3(phi_w[i-1]))*(phi_w[i-1]-phi_w[i])/2
    #int_final= int_temp + sum(intz)
    return int_temp # New value of intz: intz[i]

def f3(phi):
    return D(phi)/phi
    
#def poly20(Vwi,phi_wi): # attempt to find Q by solving a polynome (devellopment of eq.20)
#    K= rho*f*G_s(phi_wi)/( phi_b*np.pi*(R**3)*(Vwi**2) )
#    delta= np.sqrt(1+4*K*Q0)
#    res= (-1 + delta)/(2*K)
#    return K, delta, res


def PI0(phi): #(eq.7)
    return PI_ent(phi) + PI_vdw(phi) + PI_ele(phi)

def PI(phi): #(eq.7)
    if type(phi) == np.float64 or type(phi) == float :
        po=PI_ent(phi) + PI_vdw(phi) + PI_ele(phi)
        if phi>=phi_crit:
            po=PI0(phi_crit)*((1.-phi_crit/phi_cp)/(1.-phi/phi_cp))**(1./m)
    if type(phi) == np.ndarray :
        po=np.zeros(len(phi))
        for i in range (len(phi)):
            po[i]=PI_ent(phi[i]) + PI_vdw(phi[i]) + PI_ele(phi[i])
            if phi[i] >= phi_crit :
                po[i]=PI0(phi_crit)*((1.-phi_crit/phi_cp)/(1.-phi[i]/phi_cp))**(1./m)  
    return po


def PI_ent(phi): #(appendix)
    X = 6.2028*np.exp( (phi_cp - phi)*(7.9 - 3.9*(phi_cp - phi)) )
    Hall_approx = (1 + phi + phi**2 - 0.67825*phi**3 - phi**4 - 0.5*phi**5 - X*phi**6) / (1 -3*phi +3*phi**2 -1.04305*phi**3)
    return (kB*T*phi/Vp)*Hall_approx
    

def PI_vdw(phi): #(appendix)
    return (-Zn*A*phi**3)/( 48*np.pi*a**3 * (phi_cp - phi_cp**(1/3.)*phi**(2/3.))**2 )


def PI_ele(phi): #(appendix)
    epsilon= 32*np.pi*np.sqrt(2)* (3/(4.*np.pi*np.sqrt(2)))**(2/3.) # Epsilon: geometric factor
    C0= I*1000 # Initial concentration in SI units (mol/m3) (monovalent ions: I=C0)
    zeta_star = ele_charge*Z*zeta/(kB*T) # reduced zeta potential
    kappa= np.sqrt( (2*1e3*Na*(ele_charge**2)*I)/(epsr*eps0*kB*T) ) # Debye-Hückel parameter
    return (epsilon*Na*kB*T*C0/Zn)*( np.cosh( kappa*a*zeta_star/(kappa*a*phi**(-1./3.)*np.cosh(kappa*a*(1-phi**(-1./3.))) + np.sinh(kappa*a*(1-phi**(-1./3.))) ) )  -1.0 )
#------------------------------------------------------------------------------


# DATA SET USED FOR SIMULATIONS (Table 1 + addons)
# -----------------------------------------------------------------------------


# Fiber

R= 3e-2 # m (Radius R = 3 mm IN THE PAPER) 
L= 0.2 # m (Length L = 1.2 m IN THE PAPER)
Lp= 1e-9 # m/(Pa.s) (Membrane Permeability)

# Inlet condition

Um= 1.0 # m/s (Mean velocity)
Q0= Um*np.pi*(R**2) # Q0: Total flow rate for x= 0 (m3/s)
P0= 1.5e4 # Pa (Pressure P0= 100-0 kPa) 
phi_b= 1e-3 # Volumic fraction in the bulk (outside of boundary layer)

# Suspended matter

a= 1e-7 #m (Radius a = 100 nm in most cases, a= 5-1000 nm)
Vp= (4/3.)*np.pi*(a**3) #m3 (Particle volume)
zeta= 3e-2 #V (Zêta potential zeta = 30 mV)
A= 1e-20 #J (Hamaker constant)

#Suspension medium

I= 1e-5 #M (mol/l) (Ionic strength)
Z= 1.0 # Ion valence
mu_b= 1e-3 # kg/(m.s) (Dynamic viscosity)

#Physicochemical properties

T= 298.15 #K (Standard temperature)
kB= 1.381e-23 #m².kg/(s².K) (Boltzmann constant)
Na= 6.022e23 #mol-1 (Avogadro number)
rho= 1000 #kg/(m3) (Solution density of water at 298 K)
ele_charge= 1.602e-19 #C (elementary charge)
epsr= 78.4 # Water permittivity
eps0= 8.854e-12 #F/m (Vacuum permittivity)
Zn= 12.0 # Number of neighboring particles in the cell lattice
phi_cp= np.pi*np.sqrt(2)/6. # phi in an hexagonal close packing
Rm= 1./(Lp*mu_b) # m-1 (Membrane resistance)
D0 = kB*T / (6*np.pi*mu_b*a) # Simplified diffusion coefficient
m=0.3 #compressibility of the deposit
# -----------------------------------------------------------------------------


Re= rho*Um*2*R/mu_b
print("Reynolds number: Re=",Re)

if Re<= 2100:
    f= 16./Re
else:
    f= 0.0791/(Re**(1/4.))
print("Friction factor: f=",f)
print("")

phi_crit=fsolve(D_0,0.5)[0]
print ('phi_crit=',phi_crit)

phi_w=np.linspace(phi_b, phi_cp-0.1, 1000) # phi_cp is the maximal value reached by phi (close packing)

plt.figure(0)

plt.subplot(221)
plt.plot(phi_w, G_s(phi_w),label='G_s= f(phi_w)')
plt.title('Simplified function G_s (no integrations)')
plt.xlabel('Surface concentration phi_w')
plt.ylabel('G_s')
plt.legend(loc='best')

#plt.subplot(222)
#plt.plot(phi_w, G(phi_w),label='G= f(phi_w)')
#plt.title('Function G (double integration)')
#plt.xlabel('Surface concentration phi_w')
#plt.ylabel('G')
#plt.legend(loc='best')

plt.subplot(223)
plt.plot(phi_w, dPI_dphi(phi_w))
plt.title('Derivative function of the osmotic pressure PI(phi_w)')
plt.xlabel('Surface concentration phi_w')
plt.ylabel('dPI/dphi_w')
plt.legend(loc='best')

plt.subplot(224)
plt.plot(phi_w, D(phi_w))
plt.title('Diffusion coefficient D(phi_w)')
plt.xlabel('Surface concentration phi_w')
plt.ylabel('D')
plt.legend(loc='best')

plt.show()

# Plotting of PI(phi)
#------------------------------------------------------------------------------


plt.figure(1)

plt.subplot(221)
plt.plot(phi_w,PI_ent(phi_w),label='PI ent',color='C1')
plt.title('Entropic contribution to osmotic pressure: PI_ent(phi_w)')
plt.xlabel('phi_w')
plt.ylabel('PI_ent (Pa)')
plt.legend(loc='best')

plt.subplot(222)
plt.plot(phi_w,PI_vdw(phi_w),label='PI vdw',color='C2')
plt.title('Van der Walls attractive contribution to osmotic pressure: PI_vdw(phi_w)')
plt.xlabel('phi_w')
plt.ylabel('PI_vdw (Pa)')
plt.legend(loc='best')

plt.subplot(223)
plt.plot(phi_w,PI_ele(phi_w),label='PI ele',color='C3')
plt.title('Electrostatic contribution to osmotic pressure: PI_ele(phi_w)')
plt.xlabel('phi_w')
plt.ylabel('PI_ele (Pa)')
plt.legend(loc='best')

plt.subplot(224)
plt.plot(phi_w,PI(phi_w),label='PI tot')
plt.plot(phi_w,PI0(phi_w), '--',label='PI ent+vdw+elec')
plt.plot(phi_w,PI_ent(phi_w),label='PI ent')
plt.plot(phi_w,PI_vdw(phi_w),label='PI vdw')
plt.plot(phi_w,PI_ele(phi_w),label='PI ele')
plt.xlabel('phi_w')
plt.ylabel('PI (Pa)')
axes=plt.gca()
axes.set_ylim([-PI0(phi_crit),2.*PI0(phi_crit)])
plt.legend(loc='best')
plt.annotate('phi crit', xy=(phi_crit, PI(phi_crit)), xytext=(phi_crit,PI(phi_crit)/2), arrowprops=dict(facecolor='red', shrink=0.05,))


plt.title('PI total: sum off all the contributions to osmotic pressure')
plt.show()

i_max=0
for i in range(1000-1): # find phi_c as PI(phi_c)= PI maximum
    if PI(phi_w[i])>PI(phi_w[i_max]):
        i_max=i

phi_c= phi_w[i_max]
print("")
print("phi_c =",phi_c,", PI(phi_c) = PI_max =",PI(phi_c),"Pa.")
print("")
#------------------------------------------------------------------------------


# Resolution of eq. (8), (9), (20) and (22) depending on phi
#------------------------------------------------------------------------------


# Initialization (x=0)

P=[P0]
Q=[Q0]
phi_w=[phi_b]
Vw=[ (P[0] - PI(phi_w[0]))/ (Rm*mu_b) ]
dx=[0]
x=[0]
dphi= 0.01

# Loop

i=1
phi_w.append( phi_w[i-1] + dphi )
while x[i-1]<L :
#while (phi_w[i] < phi_cp - 0.1) and (x[i-1]<L):
    #print (phi_w[i])
    Vw.append( (P[i-1] - PI(phi_w[i])) / (Rm*mu_b) ) # (eq.22)
    Vwi= Vw[i]
    phi_wi= phi_w[i]
    GG= G(phi_wi) # double integration is made before the fsolve to gain computation time
    Q.append( fsolve(eq20, Q[i-1] )[0]) # (eq.20)
    dx.append( -(Q[i]- Q[i-1])/(2*np.pi*R*Vw[i]) ) # (eq.8)
    x.append( dx[i] + x[i-1])
    P.append( P[i-1] - ( dx[i]*f*rho*(Q[i])**2 / (np.pi**2 * R**5) ) )
    if dx[i]>0.03 :
        del Vw[-1]
        del Q[-1]
        del dx[-1]
        del x[-1]
        del P[-1]
        del phi_w[-1]
        dphi=dphi*0.5
        phi_w.append( phi_w[i-1] + dphi )
    else:
        i+= 1
        phi_w.append( phi_w[i-1] + dphi )
    
index=i
phi_w.remove(phi_w[index]) # value of phi_w in excess
#------------------------------------------------------------------------------


# 3D plotting of PHI=f(XX,ZZ) (XX=x and ZZ=z)
#------------------------------------------------------------------------------

#Initialisation

ZZ=[0] # X axis in the 3D plot
intz=[0] # used for ZZ
XX=[0] # Y axis in the 3D plot
PHI=[phi_b] # Z axis in the 3D plot

# Loop

for i in range(1,index):
    intz.append( inttrapz(i) )
    ZZ_res=0
    for j in range(i+1):
        XX.append( x[i] ) # XX = [x0, x1, x1, x2, x2, x2, ...]
        PHI.append( phi_w[i-j]) # PHI = [phi_b, phi_w1, phi_b, phi_w2, phi_w1, phi_b, ...]
        if j==0:
            ZZ.append(0)
        else:          
            ZZ_res=ZZ_res+ - intz[i-j+1] / Vw[i] #[0, 0, z1, 0, z1, Z2, ....]
            ZZ.append(ZZ_res)
            
            
#print(XX)
#print("")
#print(ZZ)
#print("")
#print(PHI)
#print("")
#print(intz)

x_plane, y_plane = np.meshgrid(np.linspace(0,ZZ[len(ZZ)-1],100), np.linspace(0,1.2,100))
Z = x_plane*0 + phi_crit

fig = plt.figure(2)
ax = Axes3D(fig)
ax.plot_surface(x_plane, y_plane, Z, alpha=0.2) # plot a plane at z=phi_crit
p= ax.scatter(ZZ, XX, PHI, c=PHI, cmap='hsv', linewidth=0.5)
ax.set_xlabel('z (m)')
ax.set_ylabel('x (m)')
ax.set_zlabel('phi')
plt.title('2D concentration profile: phi=f(x,z)')
fig.colorbar(p)
plt.show()
#------------------------------------------------------------------------------


print("i max=",index)
print("Size of P:",np.size(P))
print("Size of Q:",np.size(Q))
print("Size of phi_w:",np.size(phi_w))
print("Size of Vw:",np.size(Vw))
print("Size of dx:",np.size(dx))
print("Size of x:",np.size(x))
print("Size of XX:",np.size(XX))
print("Size of ZZ:",np.size(ZZ))
print("Size of PHI:",np.size(PHI))
print("Last value of x:",x[index-1]," (for phi_w=",phi_w[index-1],")")

# Plotting

plt.figure(3)

plt.subplot(221)
plt.plot(x,P,label='P(x)',color='C4')
plt.title('Transmembrane Pressure: P(x)')
plt.xlabel('x (m)')
plt.ylabel('P (Pa)')
plt.legend(loc='best')

plt.subplot(222)
plt.plot(x,Q,label='Q(x)',color='C5')
plt.title('Global flow rate: Q(x)')
plt.xlabel('x (m)')
plt.ylabel('Q (m3/s)')
plt.legend(loc='best')

plt.subplot(223)
plt.plot(x,phi_w,label='phi_w(x)',color='C6')
plt.title('Surface concentration alongside the membrane: phi_w(x)')
plt.xlabel('x (m)')
plt.ylabel('phi_w')
plt.legend(loc='best')

plt.subplot(224)
plt.plot(x,Vw,label='Vw(x)',color='C7')
plt.title('Local permeate flux: Vw(x)')
plt.xlabel('x (m)')
plt.ylabel('Vw (m/s)')
plt.legend(loc='best')

plt.show()



#plt.figure (3)
#Q20=np.linspace(0.9*Q0,Q0,100)
#plt.plot(Q20,eq20(Q20))
#plt.title("eq20=f(Q)")
#plt.show

#for i in range(0,5):
#    print("i=",i,"\n phi_w=",phi_w[i],"\n Vw=",Vw[i],"\n Q=",Q[i],"\n dx=",dx[i],"\n x=",x[i],"\n P=",P[i])
#    print("")


#Vwi=Vw[1]
#phi_wi=phi_w[1]
#print("Analytic resolution of eq.20, for Q[1]:")
#print("phi_w[1]=",phi_wi,", Vw[1]=",Vwi)
#print(" K , delta, Q:",poly20(Vwi,phi_wi))
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
finish = time.time()
print("")
print('Time spent: ', finish-start ,'s.')