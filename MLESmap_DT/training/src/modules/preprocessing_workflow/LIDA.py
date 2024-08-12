import numpy as np
from scipy.linalg import lstsq, solve
from scipy.signal import lfilter, lfiltic

def LIDA(dt,xgtt,omega,ksi,u0,ut0,rinf):
  
#%Script to Python Marisol Monterrubio 13/2/2023
#% function [u,ut,utt] = LIDA(dt,xgtt,omega,ksi,u0,ut0,rinf)
#%
#% Linear Implicit Dynamic Analysis (LIDA)
#%
#% function [u,ut,utt] = LIDA(dt,xgtt,omega,ksi,u0,ut0,rinf)
#%     Linear implicit direct time integration of second order differential
#%     equation of motion of dynamic response of linear elastic SDOF systems
#%
#% Description
#%     The General Single Step Single Solve (GSSSS) family of algorithms
#%     published by X.Zhou & K.K.Tamma (2004) is employed for direct time
#%     integration of the general linear or nonlinear structural Single
#%     Degree of Freedom (SDOF) dynamic problem. The optimal numerical
#%     dissipation and dispersion zero order displacement zero order
#%     velocity algorithm designed according to the above journal article,
#%     is used in this routine. This algorithm encompasses the scope of
#%     Linear Multi-Step (LMS) methods and is limited by the Dahlquist
#%     barrier theorem (Dahlquist,1963). The force - displacement - velocity
#%     relation of the SDOF structure is linear.
#%
#% Input parameters
#%     #dt# (scalar): time step
#%     #xgtt# ([#nstep# x 1]): column vector of the acceleration history of
#%         the excitation imposed at the base. #nstep# is the number of time
#%         steps of the dynamic response.
#%     #omega# (scalar): eigenfrequency of the structure in rad/sec
#%     #ksi# (scalar): ratio of critical damping of the SDOF system
#%     #u0# (scalar): initial displacement of the SDOF system
#%     #ut0# (scalar): initial velocity of the SDOF system
#%     #rinf# (scalar): minimum absolute value of the eigenvalues of the
#%         amplification matrix. For the amplification matrix see eq.(61) in
#%         Zhou & Tamma (2004).
#%
#% Output parameters
#%     #u# ([#nstep# x 1]): time-history of displacement
#%     #ut# ([#nstep# x 1]): time-history of velocity
#%     #utt# ([#nstep# x 1]): time-history of acceleration
#%
#% Example (Figure 6.6.1 in Chopra, Tn=1sec)
#%     dt=0.02;
#%     fid=fopen('elcentro.dat','r');
#%     text=textscan(fid,'%f %f');
#%     fclose(fid);
#%     xgtt=9.81*text{1,2};
#%     Tn=1;
#%     omega=2*pi/Tn;
#%     ksi=0.02;
#%     u0=0;
#%     ut0=0;
#%     rinf=1;
#%     [u,ut,utt] = Linear_implicit_dynamic_analysis(dt,xgtt,omega,ksi,u0,ut0,rinf);
#%     D=max(abs(u))/0.0254
#%
#%__________________________________________________________________________
#% Copyright (c) 13-Sep-2015
#%     George Papazafeiropoulos
#%     First Lieutenant, Infrastructure Engineer, Hellenic Air Force
#%     Civil Engineer, M.Sc., Ph.D. candidate, NTUA
#%     Email: gpapazafeiropoulos@yahoo.gr
#%     Website: http://users.ntua.gr/gpapazaf/
#% _________________________________________________________________________

  #%% Integration constants
  #% zero-order displacement & velocity overshooting behavior and
  #% optimal numerical dissipation and dispersion
  #print('input',dt,xgtt,omega,ksi,u0,ut0,rinf)
  w1=-15*(1-2*rinf)/(1-4*rinf); #% suggested
  w2=15*(3-4*rinf)/(1-4*rinf); #% suggested
  w3=-35*(1-rinf)/(1-4*rinf); #% suggested
  W1=(1/2+w1/3+w2/4+w3/5)/(1+w1/2+w2/3+w3/4); #% definition
  W1L1=1/(1+rinf);
  W2L2=(1/2)/(1+rinf);
  W3L3=(1/2)/(1+rinf)**2;
  W1L4=1/(1+rinf);
  W2L5=1/(1+rinf)**2; #% suggested
  W1L6=((3-rinf)/2)/(1+rinf);
  l1=1;
  l2=1/2;
  l3=(1/2)/(1+rinf);
  l4=1;
  l5=1/(1+rinf);

  #%% Transfer function denominator
  Omega=omega*dt;
  D=W1L6+2*W2L5*ksi*Omega+W3L3*Omega**2;
  A31=-Omega**2/D;
  A32=-1./D*(2.*ksi*Omega+W1L1*Omega**2);
  A33=1-1./D*(1+2.*W1L4*ksi*Omega+W2L2*Omega**2);
  A11=1+l3*A31;
  A12=l1+l3*A32;
  A13=l2-l3*(1-A33);
  A21=l5*A31;
  A22=1+l5*A32;
  A23=l4-l5*(1-A33);
  #% Amplification matrix
  A=np.asarray([[A11, A12, A13],[A21, A22, A23],[A31, A32, A33]])
  #% Amplification matrix invariants
  A1=A[0][0]+A[1][1]+A[2][2];
  A2=A[0][0]*A[1][1]-A[0][1]*A[1][0]+A[0][0]*A[2][2]-A[0][2]*A[2][0]+A[1][1]*A[2][2]-A[1][2]*A[2][1];
  A3=A[0][0]*A[1][1]*A[2][2]-A[0][0]*A[1][2]*A[2][1]-A[0][1]*A[1][0]*A[2][2]+A[0][1]*A[1][2]*A[2][0]+A[0][2]*A[1][0]*A[2][1]-A[0][2]*A[1][1]*A[2][0];
  #% Transfer function denominator
  a = np.asarray([1, -A1, A2, -A3])
  '''print('a')
  print(a)'''
  #%% Transfer function nominator
  B1=1./D*dt**2.*l3*W1;
  B2=1./D*dt**2.*(l3*(1-W1)-(A22+A33)*l3*W1+A12*l5*W1+A13*W1);
  B3=1./D*dt**2.*(-(A22+A33)*l3*(1-W1)+A12*l5*(1-W1)+A13*(1-W1)+(A22*A33-A23*A32)*l3*W1-(A12*A33-A13*A32)*l5*W1+(A12*A23-A13*A22)*W1);
  B4=1./D*dt**2.*((A22*A33-A23*A32)*l3*(1-W1)-(A12*A33-A13*A32)*l5*(1-W1)+(A12*A23-A13*A22)*(1-W1));
  
  b=np.asarray([B1,B2,B3,B4])
  '''print('b')
  print(b)'''

  #%% form initial conditions for filter function
  #% equivalent external force
  f=-xgtt;
  #print('f[0]',f,len(f),type(f))
  #% stiffness
  k=omega**2;
  #% damping constants
  c=2.*omega*ksi;
  #% initial acceleration
  utt0=-f[0]-(k*u0+c*ut0);  
  U_1=solve(A,np.asarray([u0,dt*ut0,dt**2*utt0]));
  u_1=U_1[0];
  U_2=solve(A,U_1);
  u_2=U_2[0];
  ypast=[u0,u_1,u_2];
  vinit=np.zeros((3));
  vinit = lfilter(-a[1:],np.ones(1),ypast) 

  
  #%% main dynamic analysis
  #print('types:', type(b), type(a),type(f), type(vinit))
  u=lfilter(b,a,f);
  #print('u',u[0:100],len(u))
  
  #%% calculate velocity from the following system of equations:
  #% 1st: the first scalar equation of the matrix equation (60) in [1]
  #% 2nd: equation of motion (eq.6.12.3 in Chopra:Dynamics of Structures)
  C_u=omega**2*A[0,2]*dt**2-A[0,0];  
  C_f=-A[0,2]*dt**2;
  C_ut=A[0,1]*dt-A[0,2]*dt**2*2*ksi*omega;
  
  #L=1/D*l3*dt^2*((1-W1)*[0;f(1:end-1)]+W1*f);
  L=(1/D*l3)*(dt**2)*((1-W1)*np.append(np.asarray(0),f[:-1])+W1*f);
  ut=(u+C_u*np.append(np.asarray(u0),u[0:-1])+C_f*np.append(np.asarray(0), f[0:-1])-L)/C_ut;
  #print(ut[0:100])
  #%ut=[ut0; diff(u)/dt];

  #%% calculate acceleration from equation of motion
  utt=(-omega**2)*u-2*ksi*omega*ut;
  #print(utt[0:100])
 
  return u, ut, utt

 
