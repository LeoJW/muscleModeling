% Muscle Model for testing different species and different parameters

clear all; close all;

% Force is N -> F/Fmax -> N
% Velocity is mm/s -> Lopt/s -> mm/s
% Length is mm -> L/Lopt -> mm
% k is N/mm -> dimensionless


%% Control constants

%---Primary controls

simiter = 2; % number of activation phases to compare
h = 1e-6; % step size
velBruteSize = 1e4; % number of points to solve for v
stimPhase = linspace(0.1,1,simiter); % version of tstart that varies
stimDur = linspace(0.1,0.4,simiter); % version of duration that varies
w = linspace(13,19,simiter); % cycle freq, will vary depending on species

%---Secondary controls

ncycles = 6; % number of cycles
tstart = 0.1;% point in cycle where activation begins (scaled 0 to 1)
duration = 0.1; % duration of cycle that is activated (scaled 0 to 1)

%---Hill constants

b1 = 0.25; % FLact
b2 = 0; % FLact
b = [b1,b2];
p1 = 4; % FLpas
p2 = 1; % FLpas
p = [p1,p2];

cmax = 1.8; % asymptote as v approaches -inf
vmax = 10; % maximum velocity in Lopt/s

c1 = 0.29; % from Biewener et al. (2014)
c2 = 1; % overall curvature of FV
fvc = [c1,c2,cmax,vmax];

m1 = 20; % scaling factor
m2 = 6; % horizontal translation
m3 = 0.01; % slope
m = [m1,m2,m3]; % FV curve, smooth ramp portion
delay = 50; % activation delay, in ms -> rescaled in a
gam1 = -0.98; % activation constant
gam2 = -0.94; % activation constant

%--Conversion constants

mRL = 18; %from EUST 1, TPB-EMR-dissections.xlsx, muscle resting length, mm
mtuRL = 32; % mtu resting length, mm
Lopt = mRL + 1; % mm
lamplitude = 1.2;
EMRArea = 0.0544/(0.000325*Lopt); % (mm^2) dry density in g/mm^3, mass in g
Fmax = 300e3*1e-6*EMRArea; % max force in N (convert from 300kPa to N/mm^2, multiply by EMR area)
vmaxActual = 5*Lopt; % mm/s

tslackl = mean([13.62,14.17,14.11]); % from EUST dissection on Fran's spreadsheet
tendonArea = 0.36; %(mm^2), guess based on Fran's spreadsheet
kActual = 1.6; % spring constant, N/mm, based on rat soleus loops
k = kActual*(Lopt/Fmax); % dimensionless

%---Singularity adjustments

Ftol = 0.1; % tolerance for F to avoid singularities
atol = 0.08; % tolerance for a to avoid singularities
% see FVactHinge below - added FV func to avoid singularities

for i = 1:simiter
    duration = 
    [x,v,F,wrk,pwr] = mus1(control,duration(i),w(i),C,conv,sing);
    
end
