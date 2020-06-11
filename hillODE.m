function [Fdot] = hillODE(t,F)
% ODE for MTU

x = ;
F = ;

A2 = 0.2; % amplitude of l
l = A2.*sin(w.*t) + 2; % MTU length, lmt/Lopt
ldot = A2.*w.*cos(w*t); % MTU velocity, vmt/Lopt

w = 4; % frequency in Hz or cycles/s
ncycles = 8; % number of cycles
tstart = 0.1;% point in cycle where activation begins (scaled 0 to 1)
duration = 0.4; % duration of cycle that is activated (scaled 0 to 1)

totaltime = ncycles/w; % time in s
t = linspace(0,totaltime,1e4); % time vector, 1e4 long
dt = totaltime/length(t); % time step
niter = length(t); % number of iterations in loop
lcycle = niter/ncycles; % cycle length in 1/1e4 s
startdur = round(tstart*lcycle); % start of activation in cycle
enddur = round(startdur + duration*lcycle); % duration of cycle activated in 1/1e4 s

delay = 50; % activation delay, in ms -> need to rescale in a
gam1 = -0.993; % activation constant
gam2 = -0.993; % activation constant

% Neural excitation, vector of zeros except one chunk which is 1s
ucycle = zeros(1,lcycle);
ucycle(startdur:enddur) = 1;
u = repmat(ucycle,1,ncycles);

% Activation function
d = (delay*1e-3)*(niter/totaltime); % delay, scaled
a = activationODE2(u,d,gam1,gam2);

FLactFunc = @(b,x) exp(-(((x-b(2))-1)./b(1)).^2);
b1 = 0.25; % FLact
b2 = 0; % FLact
b = [b1,b2];
FLact = FLactFunc(b,x);

p1 = 4; % FLpas
p2 = 1; % FLpas
p = [p1,p2];
FLpas = FLpasFunc(p,x)

% FV sigmoid function
Fmax = 1; % maximum force in N
cmax = 1.8; % asymptote as v approaches -inf
s1 = 1;
s2 = 1/cmax; % "asymptote", upper limit
s3 = 0.5;
s4 = 6; % affects steepness of slope at 0

% FV function for ref: FVsig = s(1)./(s(2) + z.*exp(-g*v));

xdot = (ln(((s1*Fmax.*FLact.*a)/(k(l-x)-FLpas) - s2)./s3))./s4; % dx/dt
Fdot = k(ldot - xdot); % dF/dt

end

%% More functions

%Passive force-length curve function
function [y] = FLpasFunc(p,x)
    y = p(1).*(x-p(2)).^2;
    y(x<p(2)) = 0;
end
