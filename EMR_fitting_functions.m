%% Fit functions to data in literature
FLactFunc = @(c,x) exp(-(((x-c(2))-1)./c(1)).^2);

%% Variables
hz = 500; % sample rate in samples/s
nsamp = linspace(0,1000,10e3); % number of samples
% If t is defined as linspace(0,100,10e3), work loops look normal
t = nsamp./hz; % time in seconds
w = 8.6; % cycle frequency in rad/s
A = 0.5; % amplitude of x
x = A.*sin(w.*t) + 0.5; % L/Lopt
v = A.*w.*cos(w*t); % Lengths/sec

Fmax = 18; % maximum force in N
cmax = 1.8; % asymptote as v approaches -inf
vmax = 1; % maximum velocity, normalized stretch rate?
k = 0.29; % from Biewener et al. (2014)
curv = 1; % overall curvature of FV

d = 50; % activation delay, in ms

%% Lu et al. (2011)
% Rabbit hind leg tibialis anterior
% FL Length expressed as stretch ratio
% FV using McMahon (1984) FV curve, velocity expressed as normalized
% stretch rate. Model able to predict behaviour of rabbit TA.
% Velocity expressed as normalized stretch rate
FV_lu = FV4param([k,curv,cmax,vmax],v);
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
u(300:700) = 1;
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

plot(v,Fm) % This one looks weird?

subplot(3,2,1), plot(t,u), xlabel("Time (s)"), ylabel("Neural Excitation")
subplot(3,2,2), plot(t,actvn), xlabel("Time (s)"), ylabel("Activation")
subplot(3,2,3), plot(t,x), xlabel("Time (s)"), ylabel("Length")
subplot(3,2,4), plot(t,v), xlabel("Time (s)"), ylabel("Velocity")
subplot(3,2,5), plot(x,Fm), xlabel("Length"), ylabel("Normalized Force")
subplot(3,2,6), plot(v,Fm), xlabel("Velocity"), ylabel("Normalized Force")

%% More functions

%Passive force-length curve function
function [y] = FLpasFunc(c,x)
    y = c(1).*(x-c(2)).^2;
    y(x<c(2)) = 0;
end
