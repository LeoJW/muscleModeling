%% Fit functions to data in literature
FLactFunc = @(c,x) exp(-(((x-c(2))-1)./c(1)).^2);

%% Meijer et al. (1998)
% Rat gastroc, with force normalized by maximal
% isometric force (Fce/Fceo) and length expressed as lce-lceo (i.e. peak 
% set at zero for when lce = lceo) in mm
x = linspace(-10,10,10e3);
y1 = FLactFunc([1,2],x);
plot(x,y1)
hold on;
y2 = FLactFunc([3.5,-1],x);
plot(x, y2)
hold off;
% Constants that match this data: c(1) = 3.5, c(2) = -1

%% Wakeling et al. (2012) and Lee et al. (2013) pygmy goat gastroc
% vmax = 3.59 muscle lengths/s, also used 10 muscle lengths/s in models
% cmax = 170 N
% k (curvature of FV curve) = 0.29 for fast fibres - this does not
% necessarily correspond to c(2), but was tested below.
v = linspace(-10,10,10e3);
FV_goat = FV4param([0.29,0.5,170,3.59],v);
plot(v,FV_goat)
hold on;
FV_goat_2 = FV4param([0.29,1,170,3.59],v);
plot(v,FV_goat_2)
hold off;
% Why do these look weird?

%% Lemaire et al. (2016) rat soleus
% Optimum CE length = 17.1, 21.6, 18.6 mm
% cmax (N) = 1.27, 1.25, 1.26
% a (Fmax,ce) = 0.25, 0.072, 0.20 (describes shape of FV curve)
% b (lopt,ce/s) = 0.40, 0.24, 0.40 (describes shape of FV)
% vmax = ~36 mm/s; therefore ~1.9 Lopt/s (rat 3)
% Eccentric CE force asymptote = 3.91, 3.44, 2.02
% Isometric force-slope ratio = 1.45, 0.38, 2.89
% Graph shown in paper has different x-axis

%% Josephson (1993) rat soleus
% vmax = 1.3-7 L/s, 2.8 at 30 deg C
% cmax = 162-212 kN*m^-2, 212 at 30 degrees C
% Note that rat soleus is primarily SO fibres, while EMR is FOG
% No graph to compare (FV shown is for insect flight muscle)

FV_lem = FV4param([0.5,1,1,2.8],v);
plot(v,FV_lem)
hold on;
FV_lem_2 = FV4param([0.5,1,2,2.8],v);
plot(v,FV_lem_2)
FV_lem_3 = FV4param([0.2,0.4,1.26,1.9],v);
plot(v,FV_lem_3) % This one is rat 3 of Lemaire et al. (2016). Constants a 
% and b are used for c(1) and c(2) but do not necessarily correspond.
hold off;

% Could also compare to graphs in Asmussen (1989) for rat soleus

%% Josephson (1993) rat EDL
% vmax = 13 L/s
% cmax = 209 N
% Less proportion SO, more FOG fibres than soleus
FV_EDL = FV4param([0.5,0.5,1,13],v);
plot(v,FV_EDL)
hold on;
FV_EDL_2 = FV4param([0.5,0.5,209,13],v);
plot(v,FV_EDL_2)
hold off;
% This ends up looking weird with c(3), cmax, as 1 or 209. Why?

%% Lu et al. (2011)
% Rabbit hind leg tibialis anterior
% FLpas using experimental data from Davis et al. (2003), where length is
% expressed as engineering strain and force expressed as engineering stress
FLpas_lu = FLpasFunc([4,1],x);
plot(x,FLpas_lu) % This fits curve from Chen and Zelter (1992) but it's
% based on frog muscle so not as comparable.

% FV using McMahon (1984) FV curve, velocity expressed as normalized
% stretch rate. Model able to predict behaviour of rabbit TA.
% Using k (c1) = 0.29
% Velocity expressed as normalized stretch rate
v = linspace(-10,10,10e3);
FV_lu = FV4param([0.29,1,1.8,1],v);
plot(v,FV_lu)
xlim([-1 1])
ylim([0 1.8])

%% Winters et al. (2011) rabbit TA, EDII and EDL
% EDL fast twitch, TA more fast twitch, EDII?
% Length expressed as %fibre length change
% Force can be normalized as F/Fmax
% Functions fit to curves given in paper

FLact_win_EDII = FLactFunc([0.05,0],x);
plot(x,FLact_win_EDII)
hold on;
FLact_win_EDL = FLactFunc([0.06,0],x);
plot(x,FLact_win_EDL)
FLact_win_ta = FLactFunc([0.14,0],x);
plot(x,FLact_win_ta)
hold off;

FLpas_win_EDII = FLpasFunc([40,1],x);
plot(x,FLpas_win_EDII)
hold on;
FLpas_win_EDL = FLpasFunc([30,1],x);
plot(x,FLpas_win_EDL)
FLpas_win_ta = FLpasFunc([20,1],x);
plot(x,FLpas_win_ta)
hold off;

FLpas_win_ta = FLpasFunc([20,1],x);
plot(x,FLpas_win_ta)
xlim([0.5 1.5])
ylim([0 1.2])
hold on;
FLpas_lu_ta = FLpasFunc([4,1],x);
plot(x,FLpas_lu_ta)
hold off;

FLact_win_ta = FLactFunc([0.14,0],x);
plot(x,FLact_win_ta)
xlim([0.2 2])
ylim([0 1.2])
hold on;
FLact_lu_ta = FLactFunc([0.25,0],x);
plot(x,FLact_lu_ta)
hold off;

%% Vectors for dummy work loops

%Vector for length over time
t = linspace(-10,10,10e3);
l = sin(t);
plot(t,l)

% Vector of zeros except one chunk which is 1s
u = zeros(1,10e3);
u(200:500) = 1;
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

FLtot = 18.*(FLact_lu_ta + FLpas_lu_ta); % Total FL(l)?
plot(x,FLtot)
xlim([0 2])
ylim([0 80])

% Ff = a(t)FLact(l)FV(v)
% This is the active component of muscle fibre force
% Plot against length?

Ff = blep.*FLact_lu_ta.*FV_lu;
plot(x,Ff)

% Whole muscle force is Fmax[Ff + Fp(l)]
% No need to consider pennation angle; EMR is parallel

Fm = 18.*(Ff + FLpas_lu_ta);
plot(x,Ff)

%% More functions

%Passive force-length curve function
function [y] = FLpasFunc(c,x)
    y = c(1).*(x-c(2)).^2;
    y(x<c(2)) = 0;
end
