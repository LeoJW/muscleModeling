%% EMR MTU Model

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
FV_lu = FV4param([k,curv,cmax,vmax],v);
plot(v,FV_lu)
xlim([-1 1])
ylim([0 1.8])

% Winters et al. (2011) rabbit TA
% FL passive component
FLpas_win = FLpasFunc([20,1],x);
plot(x,FLpas_win)
xlim([0.5 1.5])
ylim([0 1.2])
hold on;
FLpas_lu = FLpasFunc([4,1],x);
plot(x,FLpas_lu)
hold off;

% FL active component
FLact_win = FLactFunc([0.14,0],x);
plot(x,FLact_win)
xlim([0.2 2])
ylim([0 1.2])
hold on;
FLact_lu = FLactFunc([0.25,0],x);
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
actvn = activationODE2(u,d,-0.993,-0.993);
plot(actvn)
hold off;

%% Hill function

FLtot = Fmax.*(FLact_lu + FLpas_lu);
plot(x,FLtot)
xlim([0 2])
ylim([0 80])

% Ff = a(t)FLact(l)FV(v)
% This is the active component of muscle fibre force
% Plot against length?

Ff = actvn.*FLact_lu.*FV_lu;
plot(x,Ff)

% Whole muscle force, Fm = Fmax[Ff + Fp(l)]cos(theta) from Biewener (2014)
% No need to consider pennation angle; EMR is parallel

Fm = Fmax.*(Ff + FLpas_lu);
plot(x,Fm)
xlim([0 1])

plot(v,Fm)

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

at = [((k/M)*(l0(1)-x0(1)), zeros(1,niter-1)]; % tendon acceleration
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

% Can I approximate dFt(i) this way?

% Then need to calculate vmt(i) using vm(i) and dFt(i)

%% More functions

%Passive force-length curve function
function [y] = FLpasFunc(c,x)
    y = c(1).*(x-c(2)).^2;
    y(x<c(2)) = 0;
end
