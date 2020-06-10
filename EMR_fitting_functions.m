%% Fit functions to data in literature
clear vars;

% FL curves are normalized force plotted against L/Lopt
% FV curve is normalized force plotted against v/vmax

% FL active component
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

% Define variables

b1 = 0.25; % FLact
b2 = 0; % FLact
b = [b1,b2];
p1 = 4; % FLpas
p2 = 1; % FLpas
p = [p1,p2];
c1 = 0.29; % FV curvature of contracting phase
c2 = 1; % FV overall curvature
cmax = 1.8; % FV asymptote as x approaches -inf
vmax = 1; % FV max velocity, normalized
fvc = [c1,c2,cmax,vmax];

Fmax = 1; % max force

d = 50; % activation delay in msec
gam1 = -0.993; % activation
gam2 = -0.993; % activation
% u and a defined below

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
u = zeros(1,10e3);
u(300:800) = 1;
figure(4)
plot(u)
xlim([0 3000])
ylim([0 1.2])
hold on;

%% Activation function
gam1 = -0.993;
gam2 = -0.993;

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
plot(x,hilltest)

%% Fitting polynomials to FV and FLpas

FVpoly = polyfit(FV_lu,v,3);
FLpaspoly = polyfit(x,FLpas_lu,7);

%% Sigmoid function for FV

cmax2 = 1.8; % "asymptote" or upper limit
s = [1 1/cmax2];
vz = linspace(-5,5,1e3); % velocity
g = 6; % affects steepness of slope at 0
vo = 1; % horizontal translation
z = 0.5;

FVsig = s(1)./(s(2) + exp(-g*(vz-(1-vo))));
figure(10)
plot(vz,FVsig)
xlim([-1 1])

FVsig2 = s(1)./(s(2) + z.*exp(-g*vz));
figure(11)
plot(vz,FVsig2)
xlim([-1 1])

%% More functions

%Passive force-length curve function
function [y] = FLpasFunc(p,x)
    y = p(1).*(x-p(2)).^2;
    y(x<p(2)) = 0;
end
