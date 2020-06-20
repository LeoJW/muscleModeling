%% EMR MTU Model
clear vars;

% FL curves are normalized force plotted against L/Lopt
% velocity in (L/Lopt)/sec normalized to v/vmax
% FV curve is normalized force plotted against v/vmax

%% Hill Model

% FL active component function
FLactFunc = @(b,x) exp(-(((x-b(2))-1)./b(1)).^2);

%% Finding vm with brute force

% Variables

w = 2; % frequency in Hz or cycles/s
ncycles = 8; % number of cycles
tstart = 0.15;% point in cycle where activation begins (scaled 0 to 1)
duration = 0.6; % duration of cycle that is activated (scaled 0 to 1)

totaltime = ncycles/w; % time in s
t = linspace(0,totaltime,1e4); % time vector, 1e4 long
dt = totaltime/length(t); % time step
niter = length(t); % number of iterations in loop
lcycle = niter/ncycles; % cycle length in 1/1e4 s
startdur = round(tstart*lcycle); % start of activation in cycle
enddur = round(startdur + duration*lcycle); % duration of cycle activated in 1/1e4 s

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

delay = 50; % activation delay, in ms -> need to rescale in a
gam1 = -0.993; % activation constant
gam2 = -0.993; % activation constant

% FV curve sigmoid version
% cmax same as above
s1 = 1;
s2 = 0.5; % "asymptote", upper limit
s3 = 6; % affects steepness of slope at 0
s4 = 1/cmax;
% s = [s1,s2,s3,cmax,vmax];

% Neural excitation, vector of zeros except one chunk which is 1s
ucycle = zeros(1,lcycle);
ucycle(startdur:enddur) = 1;
u = repmat(ucycle,1,ncycles);

% Activation function
d = (delay*1e-3)*(niter/totaltime); % delay, scaled
a = activationODE2(u,d,gam1,gam2);

% Vector for Hill constants
% C = [b1,b2,p1,p2,c1,c2,cmax,vmax,Fmax];

k = 0.1; % spring constant

%% Solving for velocity with ode45

% Vector input for hillODE constants
% C2 = [b1,b2,p1,p2,s1,s2,s3,cmax,vmax,Fmax];
C = [b1,b2,p1,p2,s1,s2,s3,s4,vmax,Fmax,k];

% Additional constants
wr = 6.283185*w; % radians per second
A2 = 0.2; % amplitude of l
l = A2.*sin(wr.*t) + 2; % MTU length, l/Lopt
ldot = A2.*wr.*cos(wr*t); % MTU velocity, ldot/vmax

% Time span
tspan = [0 totaltime];

% Initial conditions
y0 = [1,k*(l(1)-1)]; % [x,F]

% Pass in an initial muscle length (x0)
% Use x0 (with l0 and ldot0) to calculate F0

% Solve ODE system
[t_out,x_out] = ode45(@(time,state) hillODE(time,t,l,ldot,state,a,C),tspan,y0);

% Boundary value solver
solinit = bvpinit(x,yinit);
sol = bvp4c(@(time,state) hillODE(time,t,l,ldot,state,a,C),bcfun,solinit);\
% bcfun is separate function specifying boundary conditions

%% More functions

%Passive force-length curve function
function [y] = FLpasFunc(p,x)
    y = p(1).*(x-p(2)).^2;
    y(x<p(2)) = 0;
end
