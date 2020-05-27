%% EMR MTU Model
clear vars;

%% Hill Model

% FL active component function
FLactFunc = @(b,x) exp(-(((x-b(2))-1)./b(1)).^2);

% Variables
hz = 500; % sample rate in samples/s
nsamp = linspace(0,1000,10e3); % number of samples
t = nsamp./hz; % time in seconds
w = 8.6; % cycle frequency in rad/s
A = 0.5; % amplitude of x
x = A.*sin(w.*t) + 0.5; % L/Lopt
v = A.*w.*cos(w*t); % Lengths/sec
plot(t,x)

b1 = 0.25; % FLact
b2 = 0; % FLact
b = [b1,b2];
p1 = 4; % FLpas
p2 = 1; % FLpas
p = [p1,p2];

Fmax = 1; % maximum force in N
cmax = 1.8; % asymptote as v approaches -inf
vmax = 10; % maximum velocity within range Wakeling (2012), Josephson (1993)
c1 = 0.29; % from Biewener et al. (2014)
c2 = 1; % overall curvature of FV
fvc = [c1,c2,cmax,vmax];

d = 50; % activation delay, in ms
gam1 = -0.993; % activation
gam2 = -0.993; % activation

% Lu et al. (2011)
FV_lu = FV4param(fvc,v); % Force-velocity
plot(v,FV_lu)
xlim([-1 1])
ylim([0 1.8])

FLpas_lu = FLpasFunc(p,x); % FL passive
plot(x,FLpas_lu)

FLact_lu = FLactFunc(b,x); % FL active
plot(x,FLact_lu)
hold off;

%% Vectors for dummy work loops

% Neural excitation, vector of zeros except one chunk which is 1s
u = zeros(1,10e3);
u(300:800) = 1;
plot(u)
xlim([0 3000])
ylim([0 1.2])
hold on;

%% Activation function
a = activationODE2(u,d,gam1,gam2);
plot(a)
hold off;

%% Hill function

% Vector for all Hill constants
C = [b1,b2,p1,p2,c1,c2,cmax,vmax,Fmax];

figure()
hilltest = hill(x,v,a,C);
plot(x,hilltest)

%% Tendon Stuff

%Constants
k = 5; % spring constant
M = 5; % mass
simTime = 10; %seconds
dt = 0.1; % time step should be small: the smaller, the more accurate
t2 = 0:dt:simTime; % time divided into time steps
niter = length(t2);

%Initial Conditions
% x0 = [2,0]; % muscle [position,velocity]
% l0 = [3,0]; % MTU [position,velocity]
% 
% %Preallocate for loop
% xm = [x0(1), zeros(1,niter-1)]; % muscle length
% vm = [x0(2), zeros(1,niter-1)]; % muscle velocity
% 
% lmt = [l0(1), zeros(1,niter-1)]; % MTU length
% vmt = [l0(2), zeros(1,niter-1)]; % MTU velocity -- we should know this as an input

% either need to solve for acceleration or pull velocity from hill model,
% or solve for velocity in some other way

% Need vm to solve for F(i) needed in minimization.
% HOW to get Ft(i) for minimization w/o solving for it first?

% Ft = [k.*(l0(1)-x0(1)), zeros(1,niter-1)]; % tendon force, equal to Fm
% dFt = [k.*(l0(2)-x0(2)), zeros(1,niter-1)]; % derivative of Ft or Fm

% dFt(i) = (k.*(lmt(i+1)-xm(i+1)-lmt(i)+xm(i)))./dt;
% Ft(i) = Ft(i-1) + dFt(i-1)*dt;
% vm(i) = vm(i-1) + am(i-1)*dt;
% vmt(i) = vmt(i-1) + amt(i-1)*dt;

%% Finding vm with minimization

% Bounds
% lb = 0;
% ub = 20;
% 
% % Initial conditions
% x0 = [2,1]; % muscle [position,velocity]
% l0 = [3,1]; % MTU [position,velocity]
% v0 = 1;

% % Preallocate
% k = 5; % spring constant
% xm = [x0(1), zeros(1,niter-1)]; % muscle length
% vm = [x0(2), zeros(1,niter-1)]; % muscle velocity
% lmt = [l0(1), zeros(1,niter-1)]; % MTU length
% vmt = [l0(2), zeros(1,niter-1)]; % MTU velocity
% Ft = [k(lmt-xm), zeros(1,niter-1)]; % tendon force, equal to Fm
% dFt = [k(vmt-vm), zeros(1,niter-1)]; % derivative of Ft or Fm

% Fvmin = Ft-hill
% Ft = kx

% Minimization

%% Finding vm with brute force

vrange = linspace(-20,20,1e4);
time = 10; % seconds
dt = 0.1; % time step
t3 = 0:dt:time;
n = length(t3);

A2 = 1; % amplitude of x
lmt = A2.*sin(w.*t) + 1; % MTU length
vmt = A2.*w.*cos(w*t); % MTU velocity

% Initial conditions
x0 = [1,1]; % muscle [xm,vm]

% Preallocate
k = 5; % spring constant
xm = [x0(1), zeros(1,n-1)]; % muscle length
vm = [x0(2), zeros(1,n-1)]; % muscle velocity
Ft = [k.*(lmt(1)-x0(1)), zeros(1,n-1)]; % tendon force
Fmin = [abs(k.*(lmt(1)-x0(1)) - hill(x0(1),x0(2),a(1),C)), zeros(1,n-1)];

% Fmin = k(l-x) - hill(x,v,a,C); want value of vrange that minimizes Fmin

for i = 2:n
    xm(i) = xm(i-1) + vm(i-1)*dt;
    Ft(i) = k.*(lmt(i-1)-xm(i-1));
    Fmin(i) = abs(k.*(lmt(i-1)-xm(i-1)) - hill(xm(i-1),vm(i-1),a(i-1),C));
end

minFmin = min(Fmin);

% find appropriate index of Fmin
% find corresponding index and value of vrange

%% More functions

%Passive force-length curve function
function [y] = FLpasFunc(p,x)
    y = p(1).*(x-p(2)).^2;
    y(x<p(2)) = 0;
end
