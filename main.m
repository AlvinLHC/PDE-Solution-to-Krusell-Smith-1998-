% Solving Krusell Smith (1998) using continuous time method and the
% extended path algorithm 
% Author: Alvin Lo Hei Chun
% Date: 20/1/2021
%==========================================================================

% close all;
clear;
clc;

tic;

%% Set Parameters 
global rho gamma alpha delta lmean 
global amin amax Na agrid da Nz zgrid ztran Zswitch amatrix zmatrix Naz

rho = 0.03;
gamma = 2;
delta = 0.025;
alpha = 1/3;
lmean = 1;

% grids 
amin = 1e-4;
amax = 100;
Na = 100;
agrid = linspace(amin,amax,Na)';
da = (amax - amin)/(Na-1);

Nz = 2;
zgrid = [0.2,1.8];
ztran = [-1/3,1/3;0.025,-0.025];
Zswitch = kron(ztran,speye(Na,Na));

amatrix = repmat(agrid,1,Nz);
zmatrix = repmat(zgrid,Na,1);
Naz = Na*Nz;

%% Solve Stationary Equilibrium 
global rSS VSS gSS KSS CSS

r0 = 0.02;
rSS = fsolve(@SolveStationaryEq,r0);
steady = stationary_eq(rSS);
VSS = steady.V;
gSS = steady.g;
KSS = steady.A;
CSS = steady.C;
YSS = KSS^(alpha)*lmean;
wSS = (1-alpha)*(alpha/(rSS+delta))^(alpha/(1-alpha));

%% Transition Dynamics 

%--------------------------------------------------------------------------
% Productivity MIT Shock
%--------------------------------------------------------------------------
global T N timespan dt Zpath
T = 200;
N = 800;
timespan = 0 + (T-0)*linspace(0,1,N)'.^1;
dt = T/N;

rhoZ = 1 - 0.95;
sigmaZ = 0.01;

Zpath = 1 + sigmaZ*exp(-rhoZ*timespan);

%--------------------------------------------------------------------------
% Conjecture Price Path 
%--------------------------------------------------------------------------
rpath = rSS*ones(N,1);

% Updating 
maxiter = 1000;
accuracy = 1e-5;
UpdateSpeed = 0.5;

TransitionDynamics;
DecomposeConsumption;

%% Plot Grath
rpath = 100*(rpath - rSS);
Kpath = 100*(log(TRAN.Kpath) - log(KSS));
wpath = 100*(log(TRAN.wpath) - log(wSS));
Cpath = 100*(log(TRAN.Cpath) - log(CSS));
Ypath = 100*(log(TRAN.Ypath) - log(YSS));
Zpath = 100*(log(Zpath) - log(1));
Paths = [rpath,Kpath,wpath,Cpath,Ypath,Zpath];
names = {'Real Interest Rate','Capital','Wage Rate','Consumption','Output','Productivity'};

figure;
for n = 1:6
    subplot(3,2,n);
    hold on;
    plot(timespan,Paths(:,n),'b','linewidth',1.5);
    plot(timespan,zeros(N,1),'k','linewidth',1);
    hold off;
    title(names(n));
    xlabel('Quarters');
    ylabel('basis point');
end

%% Plot Decomposition

C_case = 100*(log(C_case) - log(CSS));
leg = {'Base','Indirect','Direct'};

figure;
hold on;    
plot(timespan,Cpath,'o-k','linewidth',1.5);
plot(timespan,C_case(:,1),'b','linewidth',1.5);
plot(timespan,C_case(:,2),'r','linewidth',1.5);
plot(timespan,sum(C_case,2),'x-k');
plot(timespan,zeros(N,1),'k');
hold off;
legend('Base','Direct','Indirect','Direct + Indirect');
xlim([min(timespan),max(timespan)]);
ylabel('percent');
xlabel('quarters');
title('Consumption Contributions');

toc;
    
    