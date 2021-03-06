% Muscle Model (version for creating heatmaps)

clear all; close all;

% Force is N -> F/Fmax -> N
% Velocity is mm/s -> Lopt/s -> mm/s
% Length is mm -> L/Lopt -> mm
% k is N/mm -> dimensionless


%% Control constants

%---Primary controls

simiter = 10; % number of activation phases to compare
witer = 10; % number of cycle freqs to compare
h = 1e-5; % step size
velBruteSize = 1e4; % number of points to solve for v
stimPhase = linspace(0,0.78,simiter); % version of tstart that varies
w = linspace(2,20,witer); % cycle freq [avg WBF, 10 Hz for comparison]
%---Secondary controls

ncycles = 6; % number of cycles
duration = 0.2; % duration of cycle that is activated (scaled 0 to 1)

contr = [simiter,h,velBruteSize,ncycles,duration];

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

C = [b1,b2,p1,p2,cmax,vmax,c1,c2,m1,m2,m3,delay,gam1,gam2];

%--Conversion constants

mRL = 18; %from EUST 1, TPB-EMR-dissections.xlsx, muscle resting length, mm
mtuRL = 32; % mtu resting length, mm
lamplitude = 1.2;
Lopt = mRL + 0.05*mRL;
EMRArea = 0.0544/(0.000325*Lopt); % (mm^2) dry density in g/mm^3, mass in g
Fmax = 300e3*1e-6*EMRArea; % max force in N (convert from 300kPa to N/mm^2, multiply by EMR area)

tslackl = mean([13.62,14.17,14.11]); % from EUST dissection on Fran's spreadsheet
tendonArea = 0.36; %(mm^2), guess based on Fran's spreadsheet
kActual = 1.6; % spring constant, N/mm, based on rat soleus loops
k = kActual*(Lopt/Fmax); % dimensionless

conv = [mRL,mtuRL,lamplitude,EMRArea,tslackl,kActual];

%---Singularity adjustments

Ftol = 0.1; % tolerance for F to avoid singularities
atol = 0.08; % tolerance for a to avoid singularities

sing = [Ftol,atol];
% see FVactHinge below - added FV func to avoid more singularities

x = cell(1,witer);
F = cell(1,witer);
wrk = cell(1,witer);
cycNum = cell(1,witer);
simt = cell(1,witer);


%% Run simulation

% Loop through different cycle frequencies, w
tic
for i = 1:witer

    % Solve x, F and wrk for each stimPhase at each w
    [simt{i},cycNum{i},x{i},F{i},wrk{i}] = mus1(contr,stimPhase,w(i),C,conv,sing);
    
end
toc

%% 

% Make table for each cycle frequency
for i = 1:witer
    
    tbl = table(stimPhase,w,wrk);
    
end

% Define colors
hmap = heatmap(tbl,xvar,yvar,'Color Variable',wrk{1});



% figure(2)
% scatter(stimPhase,[wrk{f:simiter}],'filled')
% xlim([0 1])
% xlabel('Stimulation Phase'), ylabel('Net Work')
    
