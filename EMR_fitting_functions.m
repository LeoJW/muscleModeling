%% Fit functions to data in literature
clear vars;

FLactFunc = @(b,x) exp(-(((x-b(2))-1)./b(1)).^2);

%% Variables
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
c1 = 0.29; % from Biewener et al. (2014)
c2 = 1; % overall curvature of FV

d = 50; % activation delay, in ms

%% Lu et al. (2011)
% Rabbit hind leg tibialis anterior
% FL Length expressed as stretch ratio
% FV using McMahon (1984) FV curve, velocity expressed as normalized
% stretch rate. Model able to predict behaviour of rabbit TA.
% Velocity expressed as normalized stretch rate
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
gam1 = -0.993;
gam2 = -0.993;

a = activationODE2(u,d,gam1,gam2);
plot(a)
hold off;

%% Hill function

FLtot = Fmax.*(FLact_lu + FLpas_lu);
plot(x,FLtot)
xlim([0 2])
ylim([0 80])

% Ff = a(t)FLact(l)FV(v)
% This is the active component of muscle fibre force

Ff = a.*FLact_lu.*FV_lu;
plot(x,Ff)

% Whole muscle force, Fm = Fmax[Ff + Fp(l)]cos(theta) from Biewener (2014)
% No need to consider pennation angle; EMR is parallel

Fm = Fmax.*(Ff + FLpas_lu);
plot(x,Fm)
xlim([0 1])

plot(v,Fm)

subplot(3,2,1), plot(t,u), xlabel("Time (s)"), ylabel("Neural Excitation")
subplot(3,2,2), plot(t,a), xlabel("Time (s)"), ylabel("Activation")
subplot(3,2,3), plot(t,x), xlabel("Time (s)"), ylabel("Normalized Length")
subplot(3,2,4), plot(t,v), xlabel("Time (s)"), ylabel("Normalized Velocity")
subplot(3,2,5), plot(x,Fm), xlabel("Length"), ylabel("Normalized Force")
subplot(3,2,6), plot(v,Fm), xlabel("Velocity"), ylabel("Normalized Force")

%% Calculating work and power

% Work = Force * distance, area inside work loop
% Could be area under contraction portion minus area under lengthening
% portion, if specify period of t

wrk = trapz(x,Fm.*sign(v)); % work, area under curve w/ neg vs pos velocity
pwr = Fm.*v; % instantaneous power
subplot(3,2,1), plot(t,u), xlabel("Time (s)"), ylabel("Neural Excitation")
subplot(3,2,2), plot(t,a), xlabel("Time (s)"), ylabel("Activation")
subplot(3,2,3), plot(t,x), xlabel("Time (s)"), ylabel("Normalized Length")
subplot(3,2,4), plot(t,v), xlabel("Time (s)"), ylabel("Normalized Velocity")
subplot(3,2,5), plot(t,Fm), xlabel("Time (s)"), ylabel("Normalized Force")
subplot(3,2,6), plot(t,pwr), xlabel("Time (s)"), ylabel("Power")

%% Testing Hill function

% Define variables

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

%% More functions

%Passive force-length curve function
function [y] = FLpasFunc(p,x)
    y = p(1).*(x-p(2)).^2;
    y(x<p(2)) = 0;
end
