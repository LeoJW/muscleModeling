% KINEMATICS-DRIVEN EMR MODEL SCRIPT
clear all; close all;
tic

% Force is N -> F/Fmax -> N
% Velocity is mm/s -> Lopt/s -> mm/s
% Length is mm -> L/Lopt -> mm
% k is N/mm -> dimensionless


%**** need to fix first cycle of u/activations
%**** Need to check that timings and rates are consistent
%****X Need some method for setting up valid initial conditions
%****X Tendon length has to be defined relative to a resting length 
%**** FV4param and the handling of vmax needs to be made CONSISTENT
%**** Derivative of EMR MTU length from data needs to be continuous

%**** Constants and setup should go before loading in data: just good form
%       At least control constants should go up top. Model constants can be
%       later


%**** Time vectors I know of right now:
%       1- Kinematics data time
%       2- Not simulation time, but another time (t,l, linterp, etc. 1e4
%       long)
%       3- Simulation time (simt)

%---Primary controls *****
simiter = 3; % number of activation phases to compare
h = 1e-5; % step size (seconds!)
velBruteSize = 1e4; % number of points to solve for v
stimPhase = linspace(0.1,0.8,simiter); % version of tstart that varies

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
EMRlengthRaw = EMRa+EMRb+EMRarc; % total EMR length (mm)

% Apply Butterworth LPF, split and smooth cycles
[lTime,EMRlength,EMRcycDur,EMRfreq] = buttersplit(kineTime,EMRlengthRaw,butterOrder,butterFreq, int64(1/h));
%**** EMRlength is maybe a misleading name: probably should be MTU length

%% Constants for Hill model


%---Secondary controls

w = EMRfreq; % frequency in Hz or cycles/s
ncycles = 4; % number of cycles
tstart = 0.1;% point in cycle where activation begins (scaled 0 to 1)
duration = 0.5; % duration of cycle that is activated (scaled 0 to 1)

%---Simulation constants setup
totaltime = ncycles/w; % time in s
t = linspace(0,totaltime,1e4); % time vector %****THIS PART NEEDS REWORK
%**** dt = totaltime/length(t); % time step
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

cmax = 1.8; % asymptote as v approaches -inf
vmax = 5; % maximum velocity in Lopt/s, want to convert to mm/s

c1 = 0.29; % from Biewener et al. (2014)
c2 = 1; % overall curvature of FV
fvc = [c1,c2,cmax,vmax];

m1 = 20; % scaling factor
m2 = 6; % horizontal translation
m3 = 0.5; % slope
% m3 = 0.1;
m = [m1,m2,m3]; % FV curve, smooth ramp portion

delay = 50; % activation delay, in ms -> rescaled in a
gam1 = -0.993; % activation constant
gam2 = -0.993; % activation constant

%--Conversion constants
%Lopt = 34; % guess in mm ****** THIS HAS TO BE *MUSCLE* LENGTH, NOT MTU LENGTH
Lopt = 12.167; %*** From Bird17, WO, Morpho.xlsx 
% Fmax = 0.84; % maximum force in N  *****Calculate from stress! 300kPa*muscle area, easy
EMRArea = 0.0428/(0.000325*Lopt); %**** (mm^2) dry density in g/mm^3, mass in g
Fmax = 300e3*1e-6*EMRArea; %**** convert from 300kPa to N/mm^2, multiply by EMR area
vmaxActual = 5*Lopt; % mm/s
%**** vmax is an issue right now

% k = 0.1; % spring constant, dimensionless, want to convert to N/mm
tendonStress = 666e6; %**** (Pa, N/m^2) Random guess, anywhere from 660-1200 AFAIK
tslackl = mean([13.62,14.17,14.11])+7; %From EUST dissection on Fran's spreadsheet
tendonArea = 0.1; %(mm^2), guess based on your spreadsheet
kActual = tendonStress*1e-6*tendonArea/tslackl;
% kActual = 5; %N/mm
k = kActual*(Lopt/Fmax); % dimensionless

%---Singularity adjustments
Ftol = 0.1; % tolerance for F to avoid singularities
atol = 0.08; % tolerance for a to avoid singularities
% see FVactHinge below - added FV func to avoid singularities


%% Neural excitation and muscle activation

%---Prep figure
% close all%*** close all should really happen at the top unless you have some special reason
figure(1)
hold on
box on
grid on
col = copper(simiter); %**** col is defined multiple times!

%---Activation function
d = (delay*1e-3)*(niter/totaltime); % delay, scaled
%******OKAY. Is niter per cycle, or for the total simulation?
%            totaltime is the total simulation, but I doubt niter is

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
    plot(t,a{i},'color',col(i,:))
    xlabel("Time (s)"), ylabel("Activation")
    drawnow
    
end


%% Anonymous functions for Hill model
FLactFunc = @(b,x) exp(-(((x-b(2))-1)./b(1)).^2); % FL active component
FLpasFunc = @(p,x) heaviside(x-p(2)).*p(1).*(x-p(2)).^2; % FL passive component
FVactHinge = @(m,v) m(3)/m(1)*log(1+exp(m(1)*v-m(2))); % FV smooth ramp function


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
C = [b1,b2,p1,p2,c1,c2,cmax,vmax]; % hill

%---MTU overall length/velocity parameters
%****Go ahead and comment these top two out as well 
% wr = 2*pi*w; % frequency in radians/s
% lamplitude = 0.2; % amplitude of l
%l = lamplitude.*sin(wr.*t) + 2; % MTU length, l/Lopt
% l = 0.2*sawtooth(2*pi*t,0.2)+2; % MTU length basic asymmetric pattern

% Define length from kinematics data
EMRy = repmat(EMRlength.',1,ncycles);
% Convert time back to sec
tSec = linspace(0,max(kineTime)*EMRcycDur/length(kineTime),EMRcycDur);
EMRt = linspace(0,max(tSec)*ncycles,length(EMRy));
%**** ^What's up with tSec here? Could be done way more efficiently
%**** Also, the sample rate and duration here DON'T match the simulation


% CONVERT length and velocity to dimensionless units and prep for sim
% velocity to L/s or mm/s
linterp = interp1(EMRt,EMRy,t);
l = linterp./Lopt;
%**** HERE'S AN ISSUE. Some fixes: 
%     -Just grab fit object out of buttersplit
%     -Have only 2 sample rates in the entire script. The data sample rate,
%     and the simulation rate, defined ONLY by the step size at the very
%     top

%---Split cycles
cycL = round(length(t)/ncycles);
cycNum = repelem(1:ncycles,cycL);

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
vsweep = linspace(-1.5*vmax,1.5*vmax,velBruteSize); %*****
wrk = cell(1,simiter);
pwr = cell(1,simiter);
% Calculate FV function at all velocities
FVactVal = FV4param(fvc,vsweep);
FVhinge = FVactHinge(m,vsweep);

%---Loop through different stimulation phases
for i = 1:simiter
    
    %---Initial Conditions: Probably good to separate this out
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
%     x{i} = [1,zeros(1,length(simt)-1)]; % muscle length initial condition
    x{i} = [x0,zeros(1,length(simt)-1)];
    v{i} = [v0,zeros(1,length(simt)-1)]; % velocity
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
        % Interpolate l,a at time point *****
        tl = interp1(t,l,simt(j));
        ta = interp1(t,a{i},simt(j));
        % Solve individual components of Hill model
        FLactVal = (1-Ftol).*FLactFunc([b1,b2],x{i}(j)) + Ftol;
        % Use x(j) to solve for muscle v
        eval = k*(tl-x{i}(j)-tslackl/Lopt)*heaviside(tl-x{i}(j)-tslackl/Lopt) - ...
            (FLactVal.*(FVactVal+FVhinge).*ta + FLpasFunc([p1,p2],x{i}(j))); %****
        % Find root of function where velocity is valid
        [errval,vind] = min(abs(eval));
        err{i}(j) = eval(vind);
        v{i}(j) = vsweep(vind);
        F{i}(j) = hill(x{i}(j),v{i}(j),ta,C);
        
        % Interpolate cycle numbers
        tcycNum = interp1(t,cycNum,simt);
        % work, area under curve w/ neg vs pos velocity
        wrk{i} = -trapz(x{i}(tcycNum>3),F{i}(tcycNum>3));
        % instantaneous power
        pwr{i}(j) = F{i}(j).*v{i}(j);
        
    end
    
    %--Convert values to real units
    x{i} = x{i}*Lopt; % converts length to mm
    F{i} = F{i}*Fmax; % converts force to Newtons
    
    % Plot output
    plot(x{i}(tcycNum>2),F{i}(tcycNum>2),'color',col(i,:))
    %plot(simt,F{i},'color',col(i,:))
    drawnow
    
end

%---Aesthetics
xlabel('EMR length (mm)')
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


%****Plot Work Per Cycle
figure(4)
hold on
box on
grid on

scatter(stimPhase,[wrk{1:simiter}],'filled')
xlim([0 1])
xlabel('Stimulation Phase'), ylabel('Net Work')



%% PLOT THINGS TO HELP DEBUG


%--Muscle velocities
figure()
hold on
%MTU velocity
yyaxis right
plot(simt,gradient(interp1(t,l,simt),h),'r-')
%Muscle velocity
yyaxis left
%Muscle velocity
for i = 1:simiter
    plot(simt,v{i},'-','color',col(i,:))
end
title('Muscle Velocity')


%--Muscle length
figure()
hold on
%MTU length
yyaxis right
plot(simt,interp1(t,l,simt),'r-')
%Muscle lengths
yyaxis left
for i=1:simiter
    plot(simt,x{i}/Lopt,'-','color',col(i,:))
end
title('Muscle Length')


toc