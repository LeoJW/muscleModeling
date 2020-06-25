% KINEMATICS-DRIVEN EMR MODEL SCRIPT
clear vars; close all;

% FL curves are normalized force plotted against L/Lopt
% velocity in (L/Lopt)/sec normalized to v/vmax
% FV curve is normalized force plotted against v/vmax


%% Constants for Hill model

%---Primary controls

simiter = 6; % number of spring constants to compare
h = 1e-3; % step size
velBruteSize = 1e4; % number of points to solve for v

stimPhase = linspace(0.1,0.7,simiter); % would this be tstart?

%---Secondary controls

w = 1; % frequency in Hz or cycles/s
ncycles = 4; % number of cycles
tstart = 0.1;% point in cycle where activation begins (scaled 0 to 1)
duration = 0.4; % duration of cycle that is activated (scaled 0 to 1)

%---Simulation constants setup
totaltime = ncycles/w; % time in s
t = linspace(0,totaltime,1e4); % time vector, 1e4 long
dt = totaltime/length(t); % time step
niter = length(t); % number of iterations in loop
lcycle = niter/ncycles; % cycle length in 1/1e4 s
startdur = ceil(stimPhase*lcycle); % start of activation in cycle
enddur = ceil(startdur + duration*lcycle); % duration of cycle activated in 1/1e4 s

%---Hill constants
b1 = 0.25; % FLact
b2 = 0; % FLact
b = [b1,b2];
p1 = 4; % FLpas
p2 = 1; % FLpas
p = [p1,p2];

Fmax = 1; % maximum force in N
cmax = 1.8; % asymptote as v approaches -inf
vmax = 1; % maximum velocity within range Wakeling (2012), Josephson (1993)

s1 = 1.8; % FVsig, "asymptote" (cmax)
s4 = 1;
s2 = s1-s4; % FVsig, enforces FV(0)=1
s3 = 6; % FVsig, affects steepness of slope at 0
s = [s1,s2,s3,s4,vmax];

delay = 50; % activation delay, in ms -> rescaled in a
gam1 = -0.993; % activation constant
gam2 = -0.993; % activation constant

k = 0.1; % spring constant

%% Neural excitation and muscle activation

%---Neural excitation, vector of zeros w/ chunks of 1s
% ucycle = zeros(1,lcycle);
% ucycle(startdur:enddur) = 1;
% u = repmat(ucycle,1,ncycles);

%---Activation function
d = (delay*1e-3)*(niter/totaltime); % delay, scaled

% Prep variables for loop
ucycle = zeros(1,lcycle);
u = cell(size(stimPhase));
a = cell(size(stimPhase));

% Loop through different stimulation phases
for i = 1:simiter
    
    % Neural excitation and activation vectors
    u{i} = zeros(1,length(simt));
    a{i} = zeros(1,length(simt));
    % Loop through each time point
    for j = 1:length(simt)
        % Solve for a
        ucycle(startdur(i):enddur(i)) = 1;
        u{i}(j) = repmat(ucycle(i),1,ncycles);
        a{i}(j) = activationODE2(u{i}(j),d,gam1,gam2);
    end
    
    % Plot output
    plot(simt,a{i},'color',col(i,:))
    drawnow
    
end

%---Singularity adjustments
Ftol = 0.1; % tolerance for F to avoid singularities
atol = 0.01; % tolerance for a to avoid singularities
a = (1-atol).*a+atol;


%% Anonymous functions for Hill model

%---FL active component
FLactFunc = @(b,x) exp(-(((x-b(2))-1)./b(1)).^2);
%---FL passive component
FLpasFunc = @(p,x) heaviside(x-p(2)).*p(1).*(x-p(2)).^2;


%% Run Simulation

%---Vector input for Hill constants
B = [b1,b2,p1,p2,s1,s2,s3,s4,vmax,Fmax];

%---MTU overall length/velocity parameters
wr = 2*pi*w; % frequency in radians/s
lamplitude = 0.2; % amplitude of l
l = lamplitude.*sin(wr.*t) + 2; % MTU length, l/Lopt
ldot = lamplitude.*wr.*cos(wr*t); % MTU velocity, ldot/vmax

% Prep figure for loop
close all
figure(1)
hold on
box on
grid on
col = copper(simiter);

% Create simulation time vector
simt = 0:h:totaltime;
% Prepare variables for loop
err = cell(size(stimPhase));
Ferr = cell(size(stimPhase));
v = cell(size(stimPhase));
x = cell(size(stimPhase));
vsweep = linspace(-1,1,velBruteSize);
% Calculate FV function at all velocities
FVactVal = FVsig([s1,s2,s3,s4,vmax],vsweep);

%---Loop through different phases of activation
for i = 1:simiter
    
    % Declare vectors for simulation run
    x{i} = [1,zeros(1,length(simt)-1)]; % muscle length initial condition
    v{i} = zeros(1,length(simt)); % velocity
    err{i} = zeros(1,length(simt)); % error
    F{i} = zeros(1,length(simt)); % force
    
    % Loop through each time point
    for j = 1:length(simt)
        % Find new x from previous v
        if j~=1
            x{i}(j) = x{i}(j-1) + v{i}(j-1)*h;
        end
        % Interpolate l,a at time point
        tl = interp1(t,l,simt(j));
        ta = interp1(t,a,simt(j));
        % Solve individual components of Hill model
        FLactVal = (1-Ftol).*FLactFunc([b1,b2],x{i}(j)) + Ftol;
        % Use x(j) to solve for muscle v
        eval = k*(tl-x{i}(j)) - (FLactVal.*FVactVal.*ta + FLpasFunc([p1,p2],x{i}(j)));
        % Find root of function where velocity is valid
        [errval,vind] = min(abs(eval));
        err{i}(j) = eval(vind);
        v{i}(j) = vsweep(vind);
        
        F{i}(j) = hillv2(x{i}(j),v{i}(j),ta,B);
    end
    
    % Plot output
    plot(x{i},F{i},'color',col(i,:))
    drawnow
    
end

%---Aesthetics
xlabel('Muscle Length')
ylabel('Force')
%---Aesthetics for colorbar
colormap(copper)
cbh = colorbar;
set(cbh,'YTick',linspace(0,1,simiter))
set(cbh,'YTickLabel', num2str(k.'))

        
%% Plot error

figure(2)
hold on
box on
grid on

for i = 1:simiter
    plot(simt, err{i},'color',col(i,:))
end
%Aesthetics
title('error')
%Aesthetics for colorbar
colormap(copper)
cbh = colorbar;
set(cbh,'YTick',linspace(0,1,simiter))
set(cbh,'YTickLabel', num2str(k.'))





%% Kinematics data

kineData = readtable("2019_07_02_Tae_gut_allFlights.csv");

%---Create n*3 matrices for each point
humerus = [kineData.humerus_x, kineData.humerus_y, kineData.humerus_z];
elbow = [kineData.elbow_x, kineData.elbow_y, kineData.elbow_z];
wrist = [kineData.wrist_x, kineData.wrist_y, kineData.wrist_z];
cmc = [kineData.cmc_x, kineData.cmc_y, kineData.cmc_z];

%---Create vectors for distances between points
humerus_elbow = humerus-elbow;
elbow_wrist = elbow-wrist;
wrist_cmc = wrist-cmc;

%---Calculate angle between two vectors
ab = dot(humerus_elbow, elbow_wrist, 2);
mag_a = vecnorm(humerus_elbow,2,2);
mag_b = vecnorm(elbow_wrist,2,2);
bc = dot(elbow_wrist, wrist_cmc, 2);
mag_c = vecnorm(wrist_cmc,2,2);
% Angles in degrees
theta_elbow = acos(ab./(mag_a.*mag_b))*180/pi;
theta_wrist = acos(bc./(mag_b.*mag_c))*180/pi;


