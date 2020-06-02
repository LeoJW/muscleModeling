%% EMR MTU Model
clear vars;

% FL curves are normalized force plotted against L/Lopt
% velocity in (L/Lopt)/sec
% FV curve is normalized force plotted against v/vmax

%% Hill Model

% FL active component function
FLactFunc = @(b,x) exp(-(((x-b(2))-1)./b(1)).^2);

%% Finding vm with brute force

% Variables

time = 9.999; % seconds
dt = 0.001; % time step
t2 = 0:dt:time; % time vector, 1e4 long
w = 8; % frequency in Hz or cycles/s
niter = length(t2); % number of iterations in loop
lcycle = niter./w; % cycle length in 1/1e4 s
duration = 0.5; % duration of cycle that is activated (scaled 0 to 1)
tstart = 0.1;% point in cycle where activation begins (scaled 0 to 1)
stt2 = tstart.*lcycle;
endact = (tstart+duration).*lcycle;

b1 = 0.25; % FLact
b2 = 0; % FLact
b = [b1,b2];
p1 = 4; % FLpas
p2 = 1; % FLpas
p = [p1,p2];

Fmax = 1; % maximum force in N
cmax = 1.8; % asymptote as v approaches -inf
vmax = 1; % maximum velocity within range Wakeling (2012), Josephson (1993)
c1 = 0.29; % from Biewener et al. (2014)
c2 = 1; % overall curvature of FV
fvc = [c1,c2,cmax,vmax];

d = 50; % activation delay, in ms -> might need to rescale
gam1 = -0.993; % activation constant
gam2 = -0.993; % activation constant

% Neural excitation, vector of zeros except one chunk which is 1s
u = zeros(1,1e4);
u(stt2:endact) = 1;
u((lcycle+stt2):(lcycle+endact)) = 1;
u((2.*lcycle+stt2):(2.*lcycle+endact)) = 1;
u((3.*lcycle+stt2):(3.*lcycle+endact)) = 1;
u((4.*lcycle+stt2):(4.*lcycle+endact)) = 1;
u((5.*lcycle+stt2):(5.*lcycle+endact)) = 1;
u((6.*lcycle+stt2):(6.*lcycle+endact)) = 1;
u((7.*lcycle+stt2):(7.*lcycle+endact)) = 1;

% Activation function
a = activationODE2(u,d,gam1,gam2);

% Vector for Hill constants
C = [b1,b2,p1,p2,c1,c2,cmax,vmax,Fmax];

k = 0.1; % spring constant
vrange = linspace(-20,20,1e4); % range of possible muscle velocities

A2 = 0.2; % amplitude of lmt
lmt = A2.*sin(w.*t2) + 2; % MTU length, lmt/Lopt
vmt = A2.*w.*cos(w*t2); % MTU velocity, vmt/Lopt

% Initial conditions
x0 = [1.1461,0]; % muscle [xm,vm]

% Preallocate
xm = [x0(1), zeros(1,niter-1)]; % muscle length
vm = [x0(2), zeros(1,niter-1)]; % muscle velocity
Ft = [k.*(lmt(1)-x0(1)), zeros(1,niter-1)]; % tendon force
dFt = [k.*(vmt(1)-x0(2)), zeros(1,niter-1)]; % deriv of tendon force

% mindiff = k(l-x) - hill(x,v,a,C); want value of vrange that minimizes mindiff

for i = 2:niter
    xm(i) = xm(i-1) + vm(i-1).*dt; % problem here with vm
    Ft(i) = k.*(lmt(i)-xm(i));
    mindiff = abs(Ft(i) - hill(xm(i),vrange,a(i),C));
    [~,index] = min(mindiff); % should correspond to index of vrange that minimizes
    vm(i) = vrange(index); % value in vrange that minimizes mindiff
end

% vm(i) = vrange(index)
% find appropriate index of mindiff
% find corresponding index and value of vrange

%% More functions

%Passive force-length curve function
function [y] = FLpasFunc(p,x)
    y = p(1).*(x-p(2)).^2;
    y(x<p(2)) = 0;
end
