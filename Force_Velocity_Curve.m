clear

%% Data ([ATP] = 8, 4, 2 mM)

Data = [0		56.4048	63.6074	61.3192
    0.5		51.812	51.8626	47.4794
    1		37.4459	35.9182	31.387
    2		17.8025	13.5516	10.2112
    3		11.443	8.34	6.3895
    4		6.2643	2.8669	2.5781
    5		3.2759	1.6526	1.594
    6		2.212	1.1823	1.2117];

%% Setting up problem
N = 40; % space (strain) discretization--number of grid points in half domain
Slim = 0.18; 
dS = Slim/N;
s = (-N:1:0)*dS; % strain 
p1 = zeros(N+1,1);
p2 = zeros(N+1,1);
p3 = zeros(N+1,1);
U_NR = 0;
% State variable vector concatenates p1, p2, p2, and U_NR
PU0 = [p1; p2; p3; U_NR];

% Set metabolite concentrations, 
MgATP = [8 4 2];
MgADP = 0; 
Pi    = 0; 

% Parameters
% g = ones(1,9);
load g0
g = g0;

%% Simulating sliding and kinetics via Strang operator splitting

% moments and force
dr = 0.01; % Power-stroke Size; Units: um
kstiff1 = 1393.2; 
kstiff2 = 13275; % Non-zero velocities
% vel = -Data(:,1)'; % micron per sec
vel = -(0:0.2:6); % micron per sec

for k = 1:3

  % Zero velocity:
  [t,PU] = ode15s(@dPUdT,[0 1],PU0,[],N,dS,MgATP(k),Pi,MgADP);  
  PU = PU(end,:);
  p1 = PU(1:1*N+1);
  p2 = PU(1*N+2:2*N+2);
  p3 = PU(2*N+3:3*N+3);
  p1_0 = dS*sum(p1); p1_1 = dS*sum(s.*p1);
  p2_0 = dS*sum(p2); p2_1 = dS*sum(s.*p2);
  p3_0 = dS*sum(p3); p3_1 = dS*sum(s.*p3);
  F_active(1,k) = kstiff2*p3_0*dr + kstiff1*( p2_1 + p3_1) ;
  % Calculate XB cycling rate
  alpha3 = g(1)*283.11;
  s3 = 0.0099383;
  K_T = 0.59698; % (mM) 
  K_D = 0.194; % MgADP dissociation constant from Yamashita etal (Circ Res. 1994; 74:1027-33).
  g2 = (MgATP(k)/K_T)/(1 + MgATP(k)/K_T + MgADP/K_D);
  k3  = g2*g(5)*0.5*105.23;%;
  r_off = k3*(exp(alpha3*(s+s3).^2).*p3);
  XBCR(1,k) = dS*sum(r_off);

  % non-zero velocities
  for j = 2:length(vel)
    j
    dt = dS/abs(vel(j));
    tend = 0.5/abs(vel(j)); % ending time of simulation
    Nstep = round(tend/dt);
    % simulate kinetics for 1/2 timestep
    [t,PU] = ode15s(@dPUdT,[0 dt/2],PU0,[],N,dS,MgATP(k),Pi,MgADP);
    PU = PU(end,:); 
    for i = 1:(Nstep-1)
      % advection (sliding step)
      PU(1:1*N+0)     = PU(2:1*N+1); PU(N+1) = 0;
      PU(1*N+2:2*N+1) = PU(1*N+3:2*N+2); PU(2*N+2) = 0;
      PU(2*N+3:3*N+2) = PU(2*N+4:3*N+3); PU(3*N+3) = 0;
      % simulate kinetics for full step
      [t,PU] = ode15s(@dPUdT,[0 dt],PU,[],N,dS,MgATP(k),Pi,MgADP);
      PU = PU(end,:); 
    end
    % final advection (sliding step)
    PU(1:1*N+0)     = 0.5*(PU(2:1*N+1) + PU(1:1*N+0));         PU(N+1) = 0.5*(0 + PU(N+1));
    PU(1*N+2:2*N+1) = 0.5*(PU(1*N+3:2*N+2) + PU(1*N+2:2*N+1)); PU(2*N+2) = 0.5*(0 + PU(2*N+2));
    PU(2*N+3:3*N+2) = 0.5*(PU(2*N+4:3*N+3) + PU(2*N+3:3*N+2)); PU(3*N+3) = 0.5*(0 + PU(3*N+3));
    % final 1/2 timestep for kinetics
    [t,PU] = ode15s(@dPUdT,[0 dt/2],PU,[],N,dS,MgATP(k),Pi,MgADP);
    PU = PU(end,:);
  
    p1 = PU(1:1*N+1);
    p2 = PU(1*N+2:2*N+2);
    p3 = PU(2*N+3:3*N+3);
    p1_0 = dS*sum(p1); p1_1 = dS*sum(s.*p1);
    p2_0 = dS*sum(p2); p2_1 = dS*sum(s.*p2);
    p3_0 = dS*sum(p3); p3_1 = dS*sum(s.*p3);
    F_active(j,k) = kstiff2*p3_0*dr + kstiff1*( p2_1 + p3_1) ;
  
    % Calculate XB cycling rate
    r_off = k3*(exp(alpha3*(s+s3).^2).*p3);
    XBCR(j,k) = dS*sum(r_off);
  
  end

end

%% plots

figure(1); clf; axes('position',[0.15 0.15 0.8 0.8]); hold on;
plot(F_active(:,1), -vel,'b-','linewidth',1.5);
plot(F_active(:,2), -vel,'g-','linewidth',1.5);
plot(F_active(:,3), -vel,'r-','linewidth',1.5);
ylabel('Velocity (ML/s)','interpreter','latex','fontsize',16);
xlabel('Force (kPa)','interpreter','latex','fontsize',16);
set(gca,'fontsize',14); 
axis([0 65 0 6]);
box on;
plot(Data(:,2),Data(:,1),'bo','linewidth',1.5,'Markersize',8,'markerfacecolor',[1 1 1]);
plot(Data(:,3),Data(:,1),'go','linewidth',1.5,'Markersize',8,'markerfacecolor',[1 1 1]);
plot(Data(:,4),Data(:,1),'ro','linewidth',1.5,'Markersize',8,'markerfacecolor',[1 1 1]);
ldg = legend('8','4','2 mM');
title(ldg,'[MgATP]');

figure(2); clf; axes('position',[0.15 0.15 0.8 0.8]); hold on;
plot(-vel, F_active(:,1)'.*(-vel),'b-','linewidth',1.5);
plot(-vel, F_active(:,2)'.*(-vel),'g-','linewidth',1.5);
plot(-vel, F_active(:,3)'.*(-vel),'r-','linewidth',1.5);
xlabel('Velocity (ML/s)','interpreter','latex','fontsize',16);
ylabel('Force$\cdot$Velocity','interpreter','latex','fontsize',16);
set(gca,'fontsize',14); 
box on;
plot(Data(:,1),Data(:,2).*Data(:,1),'bo','linewidth',1.5,'Markersize',8,'markerfacecolor',[1 1 1]);
plot(Data(:,1),Data(:,3).*Data(:,1),'go','linewidth',1.5,'Markersize',8,'markerfacecolor',[1 1 1]);
plot(Data(:,1),Data(:,4).*Data(:,1),'ro','linewidth',1.5,'Markersize',8,'markerfacecolor',[1 1 1]);

figure(3); clf; axes('position',[0.15 0.15 0.8 0.8]); hold on;
plot(-vel, XBCR(:,1) ,'b-','linewidth',1.5);
plot(-vel, XBCR(:,2) ,'g-','linewidth',1.5);
plot(-vel, XBCR(:,3) ,'r-','linewidth',1.5);
xlabel('Velocity (ML/s)','interpreter','latex','fontsize',16);
ylabel('XB cycle rate (1/sec.)','interpreter','latex','fontsize',16);
set(gca,'fontsize',14); 
axis([0 6 0 25]);
box on;


