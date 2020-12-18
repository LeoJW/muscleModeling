% KINEMATICS-DRIVEN EMR MODEL SCRIPT
clear all; close all;

% Force is N -> F/Fmax -> N
% Velocity is mm/s -> Lopt/s -> mm/s
% Length is mm -> L/Lopt -> mm
% k is N/mm -> dimensionless


%% Control constants

%---Primary controls

simiter = 10; % number of activation phases to compare
h = 1e-4; % step size
stimPhase = linspace(0.1,1,simiter); % version of tstart that varies

%---Secondary controls

ncycles = 6; % number of cycles
tstart = 0.1;% point in cycle where activation begins (scaled 0 to 1)
duration = 0.5; % duration of cycle that is activated (scaled 0 to 1)

%---Hill constants

b1 = 0.25; % FLact
b2 = 0; % FLact
b = [b1,b2];
p1 = 4; % FLpas
p2 = 1; % FLpas
p = [p1,p2];

cmax = 1.8; % asymptote as v approaches -inf
vmax = 5; % maximum velocity in Lopt/s, want to convert to mm/s

c1 = 0.29; % from Biewener et al. (2014)
c2 = 1; % overall curvature of FV
fvc = [c1,c2,cmax,vmax];

m1 = 20; % scaling factor
m2 = 6; % horizontal translation
m3 = 0.01; % slope
m = [m1,m2,m3]; % FV curve, smooth ramp portion

delay = 50; % activation delay, in ms -> rescaled in a
gam1 = -0.982; % activation constant
gam2 = -0.982; % activation constant
% 0.993

%--Conversion constants

%Lopt = 12.167; % from Bird17, WO, Morpho.xlsx
Lopt = 15;
EMRArea = 0.0544/(0.000325*Lopt); % (mm^2) dry density in g/mm^3, mass in g
Fmax = 300e3*1e-6*EMRArea; % max force in N (convert from 300kPa to N/mm^2, multiply by EMR area)
vmaxActual = 5*Lopt; % mm/s

tendonE = 1000e6; % tendon elastic modulus (Pa, N/m^2), anywhere from 660-1200e6
tslackl = mean([13.62,14.17,14.11]); % from EUST dissection on Fran's spreadsheet
tendonArea = 0.36; %(mm^2), guess based on Fran's spreadsheet
kActual = tendonE*1e-6*tendonArea/tslackl; % N/mm^2
k = kActual*(Lopt/Fmax); % dimensionless

%---Singularity adjustments

Ftol = 0.1; % tolerance for F to avoid singularities
atol = 0.01; % tolerance for a to avoid singularities
% see FVactHinge below - added FV func to avoid singularities


%% Get morpho and kinematics data for EMR length

[kine,EMG] = readKineEMG();

%--Buttersplit inputs
butterOrder = 5;
butterFreq = 0.2;

% time vector for kinematics (s)
kineTime = kine(:,1);
% cycle frequency?
% Angles in deg
theta = kine(:,2); % elbow angle
phi = kine(:,3); % manus angle

% Wing geometry measurements for EMR length calculations
humL = mean([26.01,24.12,24.73]); % length of humerus
humOriginL = mean([3.17,3.81,3.66]); % how far up humerus EMR attaches, guess for now
radL = mean([32.18,31.74,32.26]); % length of radius bone
manusr = 0.5*mean([3.71,3.69,3.79]); % radius of wrist joint arc section (radius of manus)
EMRa = sqrt((radL)^2 + (humOriginL)^2 - 2*radL*humOriginL*cosd(theta));
EMRb = mean([2.37,1.98,2.02]); % how far down manus EMR attaches, fixed length
EMRarc = manusr.*(phi*pi/180); % length of EMR arc section
EMRmtuLengthRaw = EMRa+EMRb+EMRarc; % total EMR length (mm)

% Apply Butterworth LPF, split and smooth cycles
[EMRsmooth,EMRcycDur,w] = buttersplit(kineTime,EMRmtuLengthRaw,butterOrder,butterFreq);
%[processed EMR MTU length, cycle duration, cycle freq]


%% Simulation constants setup

%---Create simulation time vector
totaltime = ncycles/w; % time in s
simt = 0:h:totaltime;
modsimt = mod(length(simt),ncycles);
if modsimt>0
    simt = 0:h:totaltime+h*(ncycles-modsimt);
end
niter = length(simt); % number of iterations in loop
velBruteSize = niter; % number of points to solve for v
lcycle = round(niter/ncycles); % cycle length in 1/1e4 s
startdur = ceil(stimPhase*lcycle); % start of activation in cycle
enddur = ceil(startdur + duration*lcycle); % duration of cycle activated in 1/1e4 s


%% Neural excitation and muscle activation

%---Prep figure
figure(1)
hold on
box on
grid on
col = copper(simiter);

%---Activation function
d = round((delay*1e-3)*(niter/totaltime)); % delay, scaled

% Prep variables for loop
u = cell(niter,simiter);
a = cell(niter,simiter);

% Loop through different stimulation onsets
for i = 1:simiter
    
    % Neural excitation and activation vectors
    ucycle = zeros(1,lcycle);
    u{i} = zeros(1,niter);
    a{i} = zeros(1,niter);
    % Define each cycle
    ucycle(startdur(i):enddur(i)) = 1;
     if enddur(i)>lcycle
         ucycle(1:(enddur(i)-lcycle)) = 1;
     end
    % Define u
    u{i} = repmat(ucycle(1:lcycle),1,ncycles);
    u{i}(1:(startdur(i)-1)) = 0;
    % Solve for a
    a{i} = activationODE2(u{i},d,gam1,gam2);
    a{i} = (1-atol).*a{i}+atol;

    % Plot output
    plot(simt,a{i},'color',col(i,:))
    xlabel("Time (s)"), ylabel("Activation")
    drawnow
    
end


%% Anonymous functions for Hill model
FLactFunc = @(b,x) exp(-(((x-b(2))-1)./b(1)).^2); % FL active component
FLpasFunc = @(p,x) heaviside(x-p(2)).*p(1).*(x-p(2)).^2; % FL passive component
FVactHinge = @(m,v) m(3)/m(1)*log(1+exp(m(1)*v-m(2))); % FV smooth ramp function


%% TPB external force

TPBArea = 0.0473/(0.000325*Lopt); % (mm^2) dry density in g/mm^3, mass in g
Fmaxtpb = 300e3*1e-6*TPBArea;
tpbonset = 0.5;
tpbdur = 0.3;
starttpb = ceil(tpbonset*lcycle); % start of activation in cycle
endtpb = ceil(startdur + tpbdur*lcycle); % duration of cycle activated in 1/1e4 s
ucyctpb = zeros(1,lcycle);
ucyctpb(starttpb:endtpb) = 1;
utpb = repmat(ucyctpb,1,ncycles);
atpb = activationODE2(utpb,d,gam1,gam2);
Ftpb = Fmaxtpb*atpb;

FtpbL = 2;


%% Run Simulation

%---Vector input for Hill constants
C = [b1,b2,p1,p2,c1,c2,cmax,vmax]; % hill

%---MTU overall length/velocity parameters
tcyc = linspace(0,1,lcycle);
EMRmtuLength = EMRsmooth(tcyc); % from kinematics data
EMRy = repmat(EMRmtuLength.',1,ncycles);

% Convert length and velocity to dimensionless units and prep for sim
% velocity to L/s or mm/s
EMRmoreSmooth = fit(simt.',EMRy.','smoothingspline','SmoothingParam',0.9999999);
l = EMRmoreSmooth(simt).'./Lopt;

%---Split cycles
cycNum = repelem(1:ncycles,lcycle);

% Prep figure for loop
figure(2)
hold on
box on
grid on

% Prepare variables for loop
err = cell(1,simiter);
F = cell(1,simiter);
v = cell(1,simiter);
x = cell(1,simiter);
vsweep = linspace(-4*vmax,4*vmax,velBruteSize);
wrk = cell(1,simiter);
pwr = cell(1,simiter);
% Calculate FV function at all velocities
FVactVal = FV4param(fvc,vsweep);
FVhinge = FVactHinge(m,vsweep);

%---Loop through different stimulation phases
for i = 1:simiter
    
    % Initial Conditions
    %Velocity initial condition
    v0 = 0;
    %Find initial muscle length that is valid (assuming v0==0)
    %Sweep thru range of x0 values
    x0sweep = linspace(0,2,velBruteSize);
    %Find the one where spring and muscle forces are balanced
    x0test = k*(l(1)-x0sweep-tslackl/Lopt).*heaviside(l(1)-x0sweep-tslackl/Lopt) - ...
        (((1-Ftol).*FLactFunc([b1,b2],x0sweep)+Ftol).*a{i}(1) + FLpasFunc([p1,p2],x0sweep));
    [~,ind] = min(abs(x0test));
    x0 = x0sweep(ind);
    
    % Declare vectors for simulation run
    x{i} = [x0,zeros(1,length(simt)-1)]; % muscle length initial condition
    v{i} = [v0,zeros(1,length(simt)-1)]; % velocity initial condition
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
        % Solve individual components of Hill model
        FLactVal = (1-Ftol).*FLactFunc([b1,b2],x{i}(j)) + Ftol;
        % Use x(j) to solve for muscle v
        eval = k*(l(j)-x{i}(j)-tslackl/Lopt).*heaviside(l(j)-x{i}(j)-tslackl/Lopt) - ...
            (FLactVal.*(FVactVal+FVhinge).*a{i}(j) + FLpasFunc([p1,p2],x{i}(j)));
        % Find root of function where velocity is valid
        [errval,vind] = min(abs(eval));
        err{i}(j) = eval(vind);
        v{i}(j) = vsweep(vind);
        F{i}(j) = hill(x{i}(j),v{i}(j),a{i}(j),C);
        
        % work, area under curve w/ neg vs pos velocity
        wrk{i} = -trapz(x{i}(cycNum>(ncycles-1)),F{i}(cycNum>(ncycles-1)));
        % instantaneous power
        pwr{i}(j) = F{i}(j).*v{i}(j);
        
    end
    
    % Convert values to real units
    x{i} = x{i}*Lopt; % converts length to mm
    F{i} = F{i}*Fmax; % converts force to N
    
    % Plot output
    plot(x{i}(cycNum>(ncycles-1)),F{i}(cycNum>(ncycles-1)),'color',col(i,:))
    %plot(simt,F{i},'color',col(i,:))
    drawnow
    
end

%---Aesthetics
xlabel('EMR muscle length (mm)')
ylabel('Force (N)')
%---Aesthetics for colorbar
colormap(copper)
cbh = colorbar;
set(cbh,'YTick',linspace(0,1,simiter))
set(cbh,'YTickLabel', num2str(stimPhase.'))

        
%% Plot error

figure(3)
%***Looking at both forms of error for funzies
subplot(2,1,1)
hold on
box on
grid on
subplot(2,1,2)
hold on
box on
grid on

for i = 1:simiter
    %Instantaneous error
    subplot(2,1,1)
    plot(simt, err{i},'color',col(i,:))
    title('Instantaneous Error')
    %Cumulative error
    subplot(2,1,2)
    plot(simt, cumsum(err{i})/length(simt), 'color',col(i,:))
    title('Cumulative Error')

end
%Create colorbars
%Insant error subplot
subplot(2,1,1)
colormap(copper)
cbh = colorbar;
set(cbh,'YTick',linspace(0,1,simiter))
set(cbh,'YTickLabel', num2str(stimPhase.'))
%Cum. error subplot
subplot(2,1,2)
colormap(copper)
cbh = colorbar;
set(cbh,'YTick',linspace(0,1,simiter))
set(cbh,'YTickLabel', num2str(stimPhase.'))


%% Plot net work per cycle
figure(4)
hold on
box on
grid on

scatter(stimPhase,[wrk{1:simiter}],'filled')
xlim([0 1])
xlabel('Stimulation Phase'), ylabel('Net Work')
