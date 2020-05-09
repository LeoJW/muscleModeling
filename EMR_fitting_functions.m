%% Fit functions to data in literature
FLactFunc = @(c,x) exp(-(((x-c(2))-1)./c(1)).^2);

%% Variables
t = linspace(0,100,10e3);
w = 2; % frequency in rad/s
A = 1; % amplitude of x
x = A.*sin(w.*t) + 1;
v = A.*w.*cos(w*t);

%% Lu et al. (2011)
% Rabbit hind leg tibialis anterior
% FL Length expressed as stretch ratio
% FV using McMahon (1984) FV curve, velocity expressed as normalized
% stretch rate. Model able to predict behaviour of rabbit TA.
% Using k (c1) = 0.29 from Biewener et al. (2014)
% Velocity expressed as normalized stretch rate
FV_lu = FV4param([0.29,1,1.8,1],v);
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

% Vector of zeros except one chunk which is 1s
u = zeros(1,10e3);
u(200:250) = 1;
plot(u)
xlim([0 3000])
ylim([0 1.2])
hold on;

%% Activation function
blep = activationODE2(u,50,-0.993,-0.993);
plot(blep)
hold off;

%% Hill function
% Fmax = 18 N (Winters et al. 2011)

FLtot = 18.*(FLact_lu + FLpas_lu); % Total FL(l)?
plot(x,FLtot)
xlim([0 2])
ylim([0 80])

% Ff = a(t)FLact(l)FV(v)
% This is the active component of muscle fibre force
% Plot against length?

Ff = blep.*FLact_lu.*FV_lu;
plot(x,Ff)

% Whole muscle force is Fmax[Ff + Fp(l)]
% No need to consider pennation angle; EMR is parallel

Fm = 18.*(Ff + FLpas_lu);
plot(x,Fm)
xlim([0 1])

plot(v,Fm) % This one looks weird?

subplot(3,2,1), plot(u), xlabel("Time"), ylabel("Neural Excitation")
subplot(3,2,2), plot(blep), xlabel("Time"), ylabel("Activation")
subplot(3,2,3), plot(x), xlabel("Time"), ylabel("Length")
subplot(3,2,4), plot(v), xlabel("Time"), ylabel("Velocity")
subplot(3,2,5), plot(x,Fm), xlabel("Length"), ylabel("Normalized Force")
subplot(3,2,6), plot(v,Fm), xlabel("Velocity"), ylabel("Normalized Force")

%% More functions

%Passive force-length curve function
function [y] = FLpasFunc(c,x)
    y = c(1).*(x-c(2)).^2;
    y(x<c(2)) = 0;
end
