%% EMR MTU Model
clear vars;

%% Hill Model

% FL active component function
FLactFunc = @(c,x) exp(-(((x-c(2))-1)./c(1)).^2);

% Variables
hz = 500; % sample rate in samples/s
nsamp = linspace(0,1000,10e3); % number of samples
% If t is defined as linspace(0,100,10e3), work loops look normal
t = nsamp./hz; % time in seconds
w = 8.6; % cycle frequency in rad/s
A = 0.5; % amplitude of x
x = A.*sin(w.*t) + 0.5; % L/Lopt
v = A.*w.*cos(w*t); % Lengths/sec

Fmax = 1; % maximum force in N
cmax = 1.8; % asymptote as v approaches -inf
vmax = 10; % maximum velocity within range Wakeling (2012), Josephson (1993)
k = 0.29; % from Biewener et al. (2014)
curv = 1; % overall curvature of FV

d = 50; % activation delay, in ms

% Lu et al. (2011)
FV_lu = FV4param([k,curv,cmax,vmax],v); % Force-velocity
plot(v,FV_lu)
xlim([-1 1])
ylim([0 1.8])

FLpas_lu = FLpasFunc([4,1],x); % FL passive
plot(x,FLpas_lu)

FLact_lu = FLactFunc([0.25,0],x); % FL active
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
a = activationODE2(u,d,-0.993,-0.993);
plot(a)
hold off;

%% Hill function

% Variables

b1 = 0.25; % FLact
b2 = 0; % FLact
p1 = 4; % FLpas
p2 = 1; % FLpas
c1 = 0.29; % FV curvature of contracting phase
c2 = 1; % FV overall curvature
cmax = 1.8; % FV asymptote as x approaches -inf
vmax = 10; % FV
% u defined earlier
d = 50; % activation delay in msec
gam1 = -0.993; % activation
gam2 = -0.993; % activation
Fmax = 1;

a = activationODE2(u,d,gam1,gam2);

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
t2 = 0:dt:simTime; %<-Many ways to define time (e.g. as a linspace, or like this),
%               but code will run faster when not saving/defining things in a loop
niter = length(t2);

%Initial Conditions
x0 = [2,0]; % muscle [position,velocity]
l0 = [3,0]; % MTU [position,velocity]

%Preallocate for loop
xm = [x0(1), zeros(1,niter-1)]; % muscle length
vm = [x0(2), zeros(1,niter-1)]; % muscle velocity

lmt = [l0(1), zeros(1,niter-1)]; % MTU length
vmt = [l0(2), zeros(1,niter-1)]; % MTU velocity

at = [((k/M)*(l0(1)-x0(1))), zeros(1,niter-1)]; % tendon acceleration
amt = [diff(vmt)./dt, zeros(1,niter-1)]; % MTU acceleration
am = [k(amt-at), zeros(1,niter-1)]; % muscle acceleration

Ft = [k(lmt-xm), zeros(1,niter-1)]; % tendon force, equal to Fm
dFt = [k(vmt-vm), zeros(1,niter-1)]; % derivative of Ft or Fm
d2Ft = [diff(dFt)./dt, zeros(1,niter-1)]; % second deriv of Ft or Fm

for i = 2:niter
    xm(i) = xm(i-1) + vm(i-1)*dt;
    lmt(i) = lmt(i-1) + vmt(i-1)*dt;
    vm(i) = vm(i-1) + am(i-1)*dt;
    vmt(i) = vmt(i-1) + amt(i-1)*dt;
    Ft(i) = Ft(i-1) + dFt(i-1);
    dFt(i) = dFt(i-1) + d2Ft(i-1)*dt;
end

%% Finding vm with minimization

% Constants
niter = 100;

% Bounds
lb = 0;
ub = 20;

% Initial conditions
x0 = [2,1]; % muscle [position,velocity]
l0 = [3,1]; % MTU [position,velocity]
v0 = 1;

% Preallocate
k = 5; % spring constant
xm = [x0(1), zeros(1,niter-1)]; % muscle length
vm = [x0(2), zeros(1,niter-1)]; % muscle velocity
lmt = [l0(1), zeros(1,niter-1)]; % MTU length
vmt = [l0(2), zeros(1,niter-1)]; % MTU velocity
Ft = [k(lmt-xm), zeros(1,niter-1)]; % tendon force, equal to Fm
dFt = [k(vmt-vm), zeros(1,niter-1)]; % derivative of Ft or Fm

% Fvmin = Ft-hill
% Ft = kx

% Minimization

for i = 1:niter
    %Define function
    Fvmin = @(x,v,a,C) k*x - hill(x,v,a,C)
    %Define initial conditions
    x0
    %Minimize
    vmin = fmincon(Fvmin,v0,[],[],[],[],lb,ub);
end

%% More functions

%Passive force-length curve function
function [y] = FLpasFunc(c,x)
    y = c(1).*(x-c(2)).^2;
    y(x<c(2)) = 0;
end
