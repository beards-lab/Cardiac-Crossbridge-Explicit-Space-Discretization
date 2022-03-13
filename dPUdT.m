function f = dPUdT(~,PU,N,dS,MgATP,Pi,MgADP)
% ODE function for the d/dt operator for the cross-bridge mode.
%  first 2N-1 entries of PU represent p1(s,t)
%  second 2N-1 entries represent p2(s,t)
%  third 2N-1 entries represent p3(s,t)
%  last entry is U_NR, the fraction of myosin heads in not-relaxed state

% State Variables
p1 = PU(1:1*N+1);
p2 = PU(1*N+2:2*N+2);
p3 = PU(2*N+3:3*N+3);
U_NR = PU(3*N+4);
U_SR = 1 - U_NR;

% calculation of moments of strain distributions
s = (-N:1:0)'*dS;
p1_0 = dS*sum(p1); p1_1 = dS*sum(s.*p1);
p2_0 = dS*sum(p2); p2_1 = dS*sum(s.*p2);
p3_0 = dS*sum(p3); p3_1 = dS*sum(s.*p3);
Pu = 1 - p1_0 - p2_0 - p3_0; % unattached permissive fraction

% definition of parameters
alpha1 = 15.14;
alpha2 = 10.06;
alpha3 = 283.11;
s3 = 0.0099383;
K_Pi = 4.007;
K_T = 0.40; % (mM) 
K_D = 0.194; % MgADP dissociation constant from Yamashita etal (Circ Res. 1994; 74:1027-33).
g1 = (MgADP/K_D)/(1 + MgADP/K_D + MgATP/K_T);
g2 = (MgATP/K_T)/(1 + MgATP/K_T + MgADP/K_D);
f1 = (Pi/K_Pi)/(1 + Pi/K_Pi); f2 = 1/(1 + Pi/K_Pi); 
ka  = 330;
kd  = f1*91.02;
k1  = f2*45.491;%
k_1 = 19.395;%
k2  = 443.49;
k_2 = g1*2.9504;
k3  = g2*59.123;%;

% Force model
dr = 0.01; % Power-stroke Size; Units: um
kstiff1 = 1067.8; 
kstiff2 = 15196; 
F_active = kstiff2*p3_0*dr + kstiff1*( p2_1 + p3_1);
% F_total = F_active + 2;

% transitions between super relaxed state and non relaxed state
ksr    = 8.1352; % 
sigma0 = 30.608;
kmsr   = 250; % 

dU_NR = + ksr*(exp(F_active/sigma0))*U_SR - kmsr*U_NR*Pu  ; 
dp1   = - kd*p1 - k1*(exp(-alpha1*s).*p1) + k_1*(exp(+alpha1*s).*p2);
dp2   = + k1*(exp(-alpha1*s).*p1) - k_1*(exp(+alpha1*s).*p2) - k2*(exp(-alpha2*s).*p2) + k_2*p3;
dp3   = + k2*(exp(-alpha2*s).*p2) - k_2*p3 - k3*(exp(alpha3*(s+s3).^2).*p3);
dp1(N+1) = dp1(N+1) + ka*Pu*U_NR/dS;

f = [dp1; dp2; dp3; dU_NR];
