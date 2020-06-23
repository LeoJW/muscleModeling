%% Fit functions to data in literature
clear vars;

% FL curves are normalized force plotted against L/Lopt
% FV curve is normalized force plotted against v/vmax

% FL active component
FLactFunc = @(b,x) exp(-(((x-b(2))-1)./b(1)).^2);

%% Variables

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
c1 = 0.29; % FV curvature of contracting phase
c2 = 1; % FV overall curvature
cmax = 1.8; % FV asymptote as x approaches -inf
vmax = 2.52; % FV max velocity, normalized
fvc = [c1,c2,cmax,vmax];

Fmax = 1; % max force

delay = 50; % activation delay, in ms -> need to rescale in a
gam1 = -0.993; % activation
gam2 = -0.993; % activation
% u and a defined below

wr = 6.283185*w; % radians per second
A = 0.2; % amplitude of x
x = A.*sin(wr.*t) + 0.9; % L/Lopt
v = A.*wr.*cos(wr.*t); % Lengths/sec

%% Lu et al. (2011)
% Rabbit hind leg tibialis anterior
% FL Length expressed as stretch ratio
% FV using McMahon (1984) FV curve, velocity expressed as normalized
% stretch rate. Model able to predict behaviour of rabbit TA.
% Velocity expressed as normalized stretch rate
figure(1)
FV_lu = FV4param([c1,c2,cmax,vmax],v);
plot(v,FV_lu)
xlim([-1 1])
ylim([0 1.8])

%% Winters et al. (2011) rabbit TA
% TA fast twitch
% Length expressed as %fibre length change
% Force can be normalized as F/Fmax
% Functions fit to curves given in paper

% FL passive component
figure(2)
% FLpas_win = FLpasFunc([20,1],x);
% plot(x,FLpas_win)
% xlim([0.5 2])
% ylim([0 1.2])
FLpas_lu = FLpasFunc([4,1],x);
plot(x,FLpas_lu)

% FL active component
figure(3)
FLact_win = FLactFunc([0.14,0],x);
plot(x,FLact_win)
hold on;
FLact_lu = FLactFunc([0.25,0],x);
plot(x,FLact_lu)
hold off;

%% Vectors for dummy work loops

% Neural excitation, vector of zeros except one chunk which is 1s
ucycle = zeros(1,lcycle);
ucycle(startdur:enddur) = 1;
u = repmat(ucycle,1,ncycles);
hold on;

%% Activation function

d = (delay*1e-3)*(niter/totaltime); % delay, scaled
a = activationODE2(u,d,gam1,gam2);
figure(5)
plot(a)
hold off;

%% Hill function

figure(6)
FLtot = Fmax.*(FLact_lu + FLpas_lu);
plot(x,FLtot)

% Ff = a(t)FLact(l)FV(v)
% This is the active component of muscle fibre force

figure(7)
Ff = a.*FLact_lu.*FV_lu;
plot(x,Ff)

% Whole muscle force, Fm = Fmax[Ff + Fp(l)]cos(theta) from Biewener (2014)
% No need to consider pennation angle; EMR is parallel

figure(8)
Fm = Fmax.*(Ff + FLpas_lu);
plot(x,Fm)
xlim([0 1])

% subplot(3,2,1), plot(t,u), xlabel("Time (s)"), ylabel("Neural Excitation")
% subplot(3,2,2), plot(t,a), xlabel("Time (s)"), ylabel("Activation")
% subplot(3,2,3), plot(t,x), xlabel("Time (s)"), ylabel("Normalized Length")
% subplot(3,2,4), plot(t,v), xlabel("Time (s)"), ylabel("Normalized Velocity")
% subplot(3,2,5), plot(x,Fm), xlabel("Length"), ylabel("Normalized Force")
% subplot(3,2,6), plot(v,Fm), xlabel("Velocity"), ylabel("Normalized Force")

%% Calculating work and power

% Work = Force * distance, area inside work loop
% Could be area under contraction portion minus area under lengthening
% portion, if specify period of t

wrk = trapz(x,Fm.*sign(v)); % work, area under curve w/ neg vs pos velocity
pwr = Fm.*v; % instantaneous power
% subplot(3,2,1), plot(t,u), xlabel("Time (s)"), ylabel("Neural Excitation")
% subplot(3,2,2), plot(t,a), xlabel("Time (s)"), ylabel("Activation")
% subplot(3,2,3), plot(t,x), xlabel("Time (s)"), ylabel("Normalized Length")
% subplot(3,2,4), plot(t,v), xlabel("Time (s)"), ylabel("Normalized Velocity")
% subplot(3,2,5), plot(t,Fm), xlabel("Time (s)"), ylabel("Normalized Force")
% subplot(3,2,6), plot(t,pwr), xlabel("Time (s)"), ylabel("Power")

%% Testing Hill function

% Vector for all Hill constants
C = [b1,b2,p1,p2,c1,c2,cmax,vmax,Fmax];

figure(9)
hilltest = hill(x,v,a,C);
plot(x,hilltest), xlabel("Length"), ylabel("Force")


%% Sigmoid function for FV

s1 = 1.8; % "asymptote", upper limit (cmax)
s2 = 0.8;
s3 = 6; % affects steepness of slope at 0
s4 = 1;
% cmax and vmax same as above

s = [s1,s2,s3,s4,vmax];

FVsigtest = FVsig(s,v);
figure(10)
plot(v,FVsigtest)


%% Hill v2 with sigmoid FV

C2 = [b1,b2,p1,p2,s1,s2,s3,cmax,vmax,Fmax];
hilltest2 = hillv2(x,v,a,C2);
figure(11)
plot(x,hilltest2)

%% More functions

%Passive force-length curve function
function [y] = FLpasFunc(p,x)
    y = p(1).*(x-p(2)).^2;
    y(x<p(2)) = 0;
end
