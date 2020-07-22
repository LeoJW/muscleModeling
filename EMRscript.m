% KINEMATICS-DRIVEN EMR MODEL SCRIPT
clear vars; close all;

% Force is F/Fmax -> convert to N
% Velocity in Lopt/s -> convert to mm/sec
% Length in L/Lopt -> convert to mm
% k is dimensionless -> convert to N/mm


%% Constants for Hill model

%---Primary controls

simiter = 6; % number of activation phases to compare
h = 1e-3; % step size
velBruteSize = 1e4; % number of points to solve for v

stimPhase = linspace(0.1,0.5,simiter); % version of tstart that varies

%---Secondary controls

w = 1; % frequency in Hz or cycles/s
ncycles = 4; % number of cycles
tstart = 0.1;% point in cycle where activation begins (scaled 0 to 1)
duration = 0.5; % duration of cycle that is activated (scaled 0 to 1)

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
vmax = 1; % maximum velocity in Lopt/s, want to convert to mm/s

c1 = 0.29; % from Biewener et al. (2014)
c2 = 1; % overall curvature of FV
fvc = [c1,c2,cmax,vmax];

m1 = 20; % scaling factor
m2 = 6; % horizontal translation
m3 = 0.5; % slope
m = [m1,m2,m3]; % FV curve, smooth ramp portion

delay = 50; % activation delay, in ms -> rescaled in a
gam1 = -0.993; % activation constant
gam2 = -0.993; % activation constant

k = 0.1; % spring constant, dimensionless, want to convert to N/mm

%---Singularity adjustments
Ftol = 0.1; % tolerance for F to avoid singularities
atol = 0.08; % tolerance for a to avoid singularities
% see FVactHinge below - added FV func to avoid singularities


%% Neural excitation and muscle activation

%---Neural excitation, vector of zeros w/ chunks of 1s
% ucycle = zeros(1,lcycle);
% ucycle(startdur:enddur) = 1;
% u = repmat(ucycle,1,ncycles);

%---Prep figure
close all
figure(1)
hold on
box on
grid on
col = copper(simiter);

%---Activation function
d = (delay*1e-3)*(niter/totaltime); % delay, scaled

% Prep variables for loop
u = cell(niter,simiter);
a = cell(niter,simiter);

% Loop through different stimulation onsets
for i = 1:simiter
    
    % Neural excitation and activation vectors
    ucycle = zeros(1,lcycle);
    u{i} = zeros(1,niter);
    a{i} = zeros(1,niter);
    % Solve for a
    ucycle(startdur(i):enddur(i)) = 1;
    u{i} = repmat(ucycle,1,ncycles);
    a{i} = activationODE2(u{i},d,gam1,gam2);
    a{i} = (1-atol).*a{i}+atol;
    a{i}(niter+1:end) = [];

    % Plot output
    plot(t,a{i},'color',col(i,:))
    xlabel("Time (s)"), ylabel("Activation")
    drawnow
    
end


%---Anonymous functions for Hill model
FLactFunc = @(b,x) exp(-(((x-b(2))-1)./b(1)).^2); % FL active component
FLpasFunc = @(p,x) heaviside(x-p(2)).*p(1).*(x-p(2)).^2; % FL passive component
FVactHinge = @(m,v) m(3)/m(1)*log(1+exp(m(1)*v-m(2))); % FV smooth ramp function


%% Morpho and kinematics data

[kine,EMG] = readKineEMG();

% time vector for kinematics
kineTime = kine(:,1);
% cycle frequency?
% lengths are in mm, angles are in deg
theta = kine(:,2); % elbow angle
phi = kine(:,3); % manus angle
humActualL = [23,22,22];
humL = mean(humActualL); % length of humerus
humOriginL = (1/7)*humL; % how far up humerus EMR attaches, guess for now
radialCOR2end = [32,32,31.5];
radL = mean(radialCOR2end); % length of radius bone
manusL = 0.75*radL; % length of manus, guess for now
manusr = 0.1*manusL; % radius of wrist joint arc section (radius of manus)
% can approximate manusr with manus diameter, but guess for now
EMRa = sqrt((radL)^2 + (humOriginL)^2 - 2*radL*humOriginL*cosd(theta));
EMRb = 0.1*manusL; % how far down manus EMR attaches, fixed length, guess for now
EMRarc = manusr.*(phi*pi/180); % length of EMR arc section
EMRlength = EMRa+EMRb+EMRarc; % total EMR length
% interpolate missing values for EMR length
EMRl = fillmissing(EMRlength,'linear','SamplePoints',kineTime);



%% TPB external force

Fmaxtpb = 1;
tpbonset = 0.1;
tpbdur = 0.4;
starttpb = ceil(tpbonset*lcycle); % start of activation in cycle
endtpb = ceil(startdur + tpbdur*lcycle); % duration of cycle activated in 1/1e4 s
ucyctpb = zeros(1,lcycle);
ucyctpb(starttpb:endtpb) = 1;
utpb = repmat(ucyctpb,1,ncycles);
atpb = activationODE2(utpb,d,gam1,gam2);
Ftpb = Fmaxtpb*atpb;


%% Run Simulation

%---Vector input for Hill constants
% B = [b1,b2,p1,p2,s1,s2,s3,s4,vmax,Fmax]; % hillv2
C = [b1,b2,p1,p2,c1,c2,cmax,vmax,Fmax]; % hill

%---MTU overall length/velocity parameters
wr = 2*pi*w; % frequency in radians/s
lamplitude = 0.2; % amplitude of l
l = lamplitude.*sin(wr.*t) + 2; % MTU length, l/Lopt
ldot = lamplitude.*wr.*cos(wr*t); % MTU velocity, ldot/vmax
% can redefine with digitized kinematics data

% l = 0.4*sawtooth(2*pi*t,0.3)+1.8; % MTU length basic asymmetric pattern
% l = 2*sin(wr.*t) + 33; % example strain pattern in mm

% Prep figure for loop
figure(2)
hold on
box on
grid on
col = copper(simiter);

% Create simulation time vector
simt = 0:h:totaltime;
% Prepare variables for loop
err = cell(1,simiter);
F = cell(1,simiter);
v = cell(1,simiter);
x = cell(1,simiter);
vsweep = linspace(-vmax,vmax,velBruteSize);
wrk = cell(1,simiter);
pwr = cell(1,simiter);
% Calculate FV function at all velocities
FVactVal = FV4param(fvc,vsweep);
FVhinge = FVactHinge(m,vsweep);

%---Loop through different stimulation phases
for i = 1:simiter
    
    % Declare vectors for simulation run
    x{i} = [1,zeros(1,length(simt)-1)]; % muscle length initial condition
    v{i} = zeros(1,length(simt)); % velocity
    err{i} = zeros(1,length(simt)); % error
    F{i} = zeros(1,length(simt)); % force
    wrk{i} = zeros(1,length(simt));
    pwr{i} = zeros(1,length(simt));
    
    % Loop through each time point
    for j = 1:length(simt)
        % Find new x from previous v
        if j~=1
            x{i}(j) = x{i}(j-1) + v{i}(j-1)*h;
        end
        % Interpolate l,a at time point
        tl = interp1(t,l,simt(j));
        ta = interp1(t,a{i},simt(j));
        tldot = interp1(t,ldot,simt(j));
        % Solve individual components of Hill model
        FLactVal = (1-Ftol).*FLactFunc([b1,b2],x{i}(j)) + Ftol;
        % Use x(j) to solve for muscle v
        eval = k*(tl-x{i}(j)) - (FLactVal.*(FVactVal+FVhinge).*ta + FLpasFunc([p1,p2],x{i}(j)));
        % Find root of function where velocity is valid
        [errval,vind] = min(abs(eval));
        err{i}(j) = eval(vind);
        v{i}(j) = vsweep(vind);
        F{i}(j) = hill(x{i}(j),v{i}(j),ta,C);
        
        % Separate cycle numbers
        cycL = round(length(simt)/4);
        cycNum = [ones(1,cycL), 2*ones(1,cycL), 3*ones(1,cycL), 4*ones(1,cycL+1)];
        % work, area under curve w/ neg vs pos velocity
        wrk{i} = trapz(x{i}(cycNum>3),F{i}(cycNum>3).*-sign(v{i}(cycNum>3)));
        % instantaneous power
        pwr{i}(j) = F{i}(j).*v{i}(j);
        
    end
    
    % Separate cycle numbers
    cycL = round(length(simt)/4); 
    cycNum = [ones(1,cycL), 2*ones(1,cycL), 3*ones(1,cycL), 4*ones(1,cycL+1)];
    
    % Plot output
    plot(x{i}(cycNum>2),F{i}(cycNum>2),'color',col(i,:))
    % plot(simt,v{i},'color',col(i,:))
    drawnow
    
end

%---Aesthetics
xlabel('Muscle Length (mm)')
ylabel('Force (N)')
%---Aesthetics for colorbar
colormap(copper)
cbh = colorbar;
set(cbh,'YTick',linspace(0,1,simiter))
set(cbh,'YTickLabel', num2str(stimPhase.'))

        
%% Plot error

figure(3)
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
set(cbh,'YTickLabel', num2str(stimPhase.'))

figure(4)
hold on
box on
grid on

bob = [wrk{1:6}];
scatter(stimPhase,bob,'filled')
xlim([0 1])
xlabel('Stimulation Phase'), ylabel('Net Work')

figure(5)
hold on
box on
grid on
plot(kineTime,EMRl)
xlabel('Time (s)'), ylabel('EMR MTU Length (mm)')


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


