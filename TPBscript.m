% KINEMATICS-DRIVEN EMR MODEL SCRIPT - WITH TPB ADDED
clear all; close all;

% Force is N -> F/Fmax -> N
% Velocity is mm/s -> Lopt/s -> mm/s
% Length is mm -> L/Lopt -> mm
% k is N/mm -> dimensionless


%% Control constants

%---Primary controls

simiter = 2; % number of activation phases to compare
h = 0.5e-5; % step size
stimPhase = linspace(0.1,1,simiter); % version of tstart that varies
velBruteSize = 1e5; % number of points to evaluate for root finding

%---Secondary controls

ncycles = 6; % number of cycles
tstart = 0.1;% point in cycle where activation begins (scaled 0 to 1)
duration = 0.4; % duration of cycle that is activated (scaled 0 to 1)

%---Hill constants

b1 = 0.25; % FLact
b2 = 0; % FLact
b = [b1,b2];
p1 = 4; % FLpas
p2 = 1; % FLpas
p = [p1,p2];

cmax = 1.8; % asymptote as v approaches -inf
% vmax = 5; % maximum velocity in Lopt/s
vmax = 5*4; % maximum velocity in Lopt/s

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
Lopt = mean([18.48,17.88,18.57]); % from Fran spreadsheet EUST1
lopt1 = mean([5.33,4.89,5.85]); % from spreadsheet EUST1, based on TPB insertion
lopt2 = Lopt - lopt1; % from spreadsheet EUST1, based on TPB insertion
EMRArea = 0.0544/(0.000325*Lopt); % (mm^2) dry density in g/mm^3, mass in g
Fmax = 300e3*1e-6*EMRArea; % max force in N (convert from 300kPa to N/mm^2, multiply by EMR area)
vmaxActual = 5*Lopt; % mm/s

tendonE = 1e9; % tendon elastic modulus (Pa, N/m^2), anywhere from 60.66-1.2e9
tslackl = mean([13.62,14.17,14.11]); % from EUST dissection on Fran's spreadsheet
tendonArea = 0.36; %(mm^2), guess based on Fran's spreadsheet
kActual = tendonE*1e-6*tendonArea/tslackl; % N/mm^2
k = kActual*(1/Fmax); % (1/mm)

%---Singularity adjustments

Ftol = 0.1; % tolerance for F to avoid singularities
atol = 0.08; % tolerance for a to avoid singularities
precision = 20; %significant digits to round to for floating-point errors
% see FVactHinge below - added FV func to avoid singularities


%% Get morpho and kinematics data for EMR straight MTU length

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
%[EMRsmooth,EMRcycDur,w] = buttersplit(kineTime,EMRmtuLengthRaw,butterOrder,butterFreq);
%[processed EMR MTU length, cycle duration, cycle freq]


%% Simulation constants setup

%---Create simulation time vector
w = 18; % cycle freq in Hz
totaltime = ncycles/w; % time in s
simt = 0:h:totaltime;
modsimt = mod(length(simt),ncycles);
if modsimt>0
    simt = 0:h:totaltime+h*(ncycles-modsimt);
end
niter = length(simt); % number of iterations in loop
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
    a{i} = activationODE2(u{i},d,gam1,gam2,1/h);
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

Lopttpb = mean([21.95,22.04,21.7]);
TPBArea = 0.0473/(0.000325*Lopttpb); % (mm^2) dry density in g/mm^3, mass in g 
Fmaxtpb = 300e3*1e-6*TPBArea;
tpbonset = 0.5;
tpbdur = 0.3;
starttpb = ceil(tpbonset*lcycle); % start of activation in cycle
endtpb = ceil(starttpb + tpbdur*lcycle); % duration of cycle activated in 1/1e4 s
ucyctpb = zeros(1,lcycle);
ucyctpb(starttpb:endtpb) = 1;
utpb = repmat(ucyctpb,1,ncycles);
atpb = activationODE2(utpb,d,gam1,gam2,1/h);
Ftpb = Fmaxtpb*atpb/Fmax;

Ftpb = zeros(size(Ftpb));


%% Run Simulation

%---Vector input for Hill constants
C = [b1,b2,p1,p2,c1,c2,cmax,vmax]; % hill

%---MTU overall length/velocity parameters
% tcyc = linspace(0,1,lcycle);
% EMRmtuLength = EMRsmooth(tcyc); % from kinematics data
% EMRy = repmat(EMRmtuLength.',1,ncycles);

% Convert length and velocity to dimensionless units and prep for sim
% velocity to L/s or mm/s
% EMRmoreSmooth = fit(simt.',EMRy.','smoothingspline','SmoothingParam',0.9999999);
%L = EMRmoreSmooth(simt).'; %<- note: No ./ needed, just / is fine

% what happens if using sin pattern instead
w = 18; % cycle frequency in Hz
wr = 2*pi*w; % convert to radians
Lamplitude = 2;
L = Lamplitude.*sin(wr.*simt) + 35; % MTU length in mm
Ldot = Lamplitude*wr.*sin(wr.*simt); % MTU velocity in mm/s

%---Split cycles
cycNum = repelem(1:ncycles,lcycle);

% Prep figure for loop
figure(2)
hold on
box on
grid on

% Prepare variables for loop to solve for l2
err2 = cell(1,simiter);
err1 = cell(1,simiter);
F2 = cell(1,simiter);
F1 = cell(1,simiter);
v2 = cell(1,simiter);
v1 = cell(1,simiter);
l2 = cell(1,simiter);
l1 = cell(1,simiter);
angle2 = cell(1,simiter);
angle1 = cell(1,simiter);
angle1v = cell(1,simiter);
angle2v = cell(1,simiter);
lt = cell(1,simiter);
ltdot = cell(1,simiter);
% vsweep = linspace(-4*vmax,4*vmax,velBruteSize);
vsweep = linspace(-vmax,vmax,velBruteSize);
wrk = cell(1,simiter);
pwr = cell(1,simiter);
% Calculate FV function at all velocities
FVactVal = FV4param(fvc,vsweep);
FVhinge = FVactHinge(m,vsweep);

%---Loop through different stimulation phases
tic
parfor i = 1:simiter
    
    % Initial Conditions
    
    %Find initial muscle length that is valid (assuming v0==0)
    %Sweep thru range of x0 values
    x0sweep = linspace(0,2,velBruteSize)*Lopt;
    %Find the one where spring and muscle forces are balanced
    x0test = k*(L(1)-x0sweep-tslackl).*heaviside(L(1)-x0sweep-tslackl) - ...
        (((1-Ftol).*FLactFunc([b1,b2],x0sweep/Lopt)+Ftol).*a{i}(1) + FLpasFunc([p1,p2],x0sweep/Lopt));
    [~,ind] = min(abs(x0test)); % here k is in 1/mm so everything becomes normalized
    x0 = x0sweep(ind)/Lopt; % nondimensional to apply to lopt1 and lopt2
    
    
    %Velocity initial condition
    v20 = 0;
    v10 = 0;
    ltdot0 = 0;
    angle1v0 = 0;
    angle2v0 = 0;
    %Angle initial conditions
    angle10 = 0;
    angle20 = 0;
    %Set initial muscle and tendon length
    l10 = lopt1*x0; % initial length of muscle section 1 at rest, mm
    l20 = lopt2*x0; % initial length of muscle section 2 at rest, mm
    lt0 = L(1) - l10 - l20; % initial tendon length, mm
    %Get intial force (same for all elements as initial angles are 0)
    F0 = k*(L(1)-x0*Lopt-tslackl).*heaviside(L(1)-x0*Lopt-tslackl); % dimensionless, F/Fmax
    
    
    % Declare vectors for simulation run
    l1{i} = [l10, zeros(1,length(simt)-1)]; % muscle length section 1, mm
    l2{i} = [l20, zeros(1,length(simt)-1)]; % muscle length section 2, mm
    lt{i} = [lt0, zeros(1,length(simt)-1)]; % tendon length, mm
    angle1{i} = [angle10, zeros(1,length(simt)-1)]; % angle 1
    angle2{i} = [angle20, zeros(1,length(simt)-1)]; % angle 2
    v1{i} = [v10, zeros(1,length(simt)-1)]; % velocity section 1
    v2{i} = [v20, zeros(1,length(simt)-1)]; % velocity section 2
    ltdot{i} = [ltdot0, zeros(1,length(simt)-1)]; % velocity of tendon
    angle1v{i} = [angle1v0, zeros(1,length(simt)-1)]; % velocity angle 1
    angle2v{i} = [angle2v0, zeros(1,length(simt)-1)]; % velocity angle 1
    F2{i} = [F0, zeros(1,length(simt)-1)]; % force muscle section 2
    F1{i} = [F0, zeros(1,length(simt)-1)]; % force muscle section 1
    err2{i} = zeros(1,length(simt)); % error
    err1{i} = zeros(1,length(simt)); % error
    wrk{i} = zeros(1,length(simt));
    pwr{i} = zeros(1,length(simt));
    
    % Loop through each time point
    for j = 2:length(simt)
        
        % Define state variables at time step (j) from previous derivatives (j-1)
        % angle1, angle2
        angle1{i}(j) = angle1{i}(j-1) + angle1v{i}(j-1)*h;
        angle2{i}(j) = angle2{i}(j-1) + angle2v{i}(j-1)*h;
        % l1, l2 
        l2{i}(j) = l2{i}(j-1) + v2{i}(j-1)*h; % mm
        l1{i}(j) = l1{i}(j-1) + v1{i}(j-1)*h; % mm
        % lt
        lt{i}(j) = lt{i}(j-1) + ltdot{i}(j-1)*h;
        
        
        % Now use current values to define derivatives at t=j
        

        % Solve individual components of Hill model
        FLactVal2 = (1-Ftol).*FLactFunc([b1,b2],l2{i}(j)/lopt2) + Ftol;
        FLactVal1 = (1-Ftol).*FLactFunc([b1,b2],l1{i}(j)/lopt1) + Ftol;
        % Use l2(j) to solve for muscle section 2 v
        eval = k*(lt{i}(j)-tslackl).*heaviside(lt{i}(j)-tslackl) - ...
            (FLactVal2.*(FVactVal+FVhinge).*a{i}(j) + FLpasFunc([p1,p2],l2{i}(j)/lopt2)); % normalized, F/Fmax
        % Find root of function where velocity is valid
        [errval,v2ind] = min(abs(eval));
        err2{i}(j) = eval(v2ind);
        v2{i}(j) = vsweep(v2ind)*lopt2; %mm/s
        F2{i}(j) = hill(l2{i}(j)/lopt2, v2{i}(j)/lopt2, a{i}(j), C); %F/Fmax
        
        %Solve for muscle section 1 v using Ftpb equation
        if Ftpb(j) < 0.01
            % Set v1 just from v2 when no Y forces
            v1{i}(j) = v2{i}(j)*lopt1/lopt2; %mm/s
            F1{i}(j) = hill(l1{i}(j)/lopt1, v1{i}(j)/lopt1, a{i}(j), C); % F/Fmax
        else
            evalagain = Ftpb(j) - F2{i}(j).*cos(angle2{i}(j)).*tan(angle1{i}(j)) - ...
                (FLactVal1.*(FVactVal+FVhinge).*a{i}(j) + ...
                FLpasFunc([p1,p2],l1{i}(j)/lopt1)).*cos(angle1{i}(j)).*tan(angle2{i}(j)); % F/Fmax
            % Find root of function where velocity is valid
            [errvalagain,v1ind] = min(abs(evalagain));
            err1{i}(j) = evalagain(v1ind);
            % Use root to find v1
            v1{i}(j) = vsweep(v1ind)*lopt1; % mm/s
            % Use v1 to find F1
            F1{i}(j) = hill(l1{i}(j)/lopt1, v1{i}(j)/lopt1, a{i}(j), C); % F/Fmax
        end
        
        
        % Calculate lt and ltdot from l1, l2, v1, v2
        ltdot{i}(j) = (Ldot(j) - v1{i}(j)*cos(angle1{i}(j)) - v2{i}(j)*cos(angle2{i}(j)))/cos(angle2{i}(j)); % mm/s
        
        
        angle1v{i}(j) = -( ...
            ( ((2*(l2{i}(j)+lt{i}(j))*(v2{i}(j) + ltdot{i}(j)) + 2*L(j)*Ldot(j) - 2*l1{i}(j)*v1{i}(j)))/(2*L(j)*(l2{i}(j)+lt{i}(j))) ...
            - ((v2{i}(j) + ltdot{i}(j))*( (l2{i}(j)+lt{i}(j))^2 + L(j)^2 - l1{i}(j)^2 ))/(2*L(j)*(l2{i}(j)+lt{i}(j))^2) ...
            - (Ldot(j)*((l2{i}(j)+lt{i}(j))^2 + L(j)^2 - l1{i}(j)^2))/(2*L(j)^2*(l2{i}(j)+lt{i}(j))) ) / ...
            (sqrt(1 - ( ((l2{i}(j)+lt{i}(j))^2 + L(j)^2 -l1{i}(j)^2)^2/(4*L(j)^2*(l2{i}(j)+lt{i}(j))^2) )) ) ...
            );
        
        angle2v{i}(j) = -( ...
            ((-2*(l2{i}(j)+lt{i}(j))*(v2{i}(j)+ltdot{i}(j)) + 2*L(j)*Ldot(j) + 2*l1{i}(j)*v1{i}(j)) / (2*L(j)*l1{i}(j)) ...
            - ( (Ldot(j)*(-(l2{i}(j)+lt{i}(j))^2 + L(j)^2 + l1{i}(j)^2))/(2*L(j)^2*l1{i}(j)) ) ...
            - ( (v1{i}(j)*(-(l2{i}(j)+lt{i}(j))^2 + L(j)^2 + l1{i}(j)^2))/(2*L(j)*l1{i}(j)^2) )) / ...
            (sqrt(1 - (( -(l2{i}(j)+lt{i}(j))^2 + L(j)^2 + l1{i}(j)^2 )^2 / (4*L(j)^2*l1{i}(j)^2) ))) ...
            );
        
        if Ftpb(j) < 0.01
            angle1v{i}(j) = 0;
            angle2v{i}(j) = 0;
        end
        
        % Calculate instantaneous power
        pwr{i}(j) = F2{i}(j).*v2{i}(j);
    end

    % work, area under curve w/ neg vs pos velocity
    % will need to specify which sections of muscle we are solving for
    wrk{i} = -trapz(l2{i}(cycNum>(ncycles-1)),F2{i}(cycNum>(ncycles-1))); % should these be calculated 
    
    % Convert values to real units
    F2{i} = F2{i}*Fmax; % converts force to N
    F1{i} = F1{i}*Fmax; % converts force to N
    
end
toc

% Loop again to plot
for i = 1:simiter
    % Plot output
    plot(l2{i}(cycNum>(ncycles-1)),F2{i}(cycNum>(ncycles-1)),'color',col(i,:))
    %plot(simt,F{i},'color',col(i,:))
end
%---Aesthetics
xlabel('EMR muscle length (mm)')
ylabel('Force (N)')
%---Aesthetics for colorbar
colormap(copper)
cbh = colorbar;
set(cbh,'YTick',linspace(0,1,simiter))
set(cbh,'YTickLabel', num2str(stimPhase.'))

        
% %% Plot error
% 
% figure(3)
% %***Looking at both forms of error for funzies
% subplot(2,1,1)
% hold on
% box on
% grid on
% subplot(2,1,2)
% hold on
% box on
% grid on
% 
% for i = 1:simiter
%     %Instantaneous error
%     subplot(2,1,1)
%     plot(simt, err2{i},'color',col(i,:))
%     title('Instantaneous Error')
%     %Cumulative error
%     subplot(2,1,2)
%     plot(simt, cumsum(err2{i})/length(simt), 'color',col(i,:))
%     title('Cumulative Error')
% 
% end
% %Create colorbars
% %Insant error subplot
% subplot(2,1,1)
% colormap(copper)
% cbh = colorbar;
% set(cbh,'YTick',linspace(0,1,simiter))
% set(cbh,'YTickLabel', num2str(stimPhase.'))
% %Cum. error subplot
% subplot(2,1,2)
% colormap(copper)
% cbh = colorbar;
% set(cbh,'YTick',linspace(0,1,simiter))
% set(cbh,'YTickLabel', num2str(stimPhase.'))
% 
% 
% %% Plot net work per cycle
% figure(4)
% hold on
% box on
% grid on
% 
% scatter(stimPhase,[wrk{1:simiter}],'filled')
% xlim([0 1])
% xlabel('Stimulation Phase'), ylabel('Net Work')


%% Look at some stuff


figure()
subplot(3,1,1)
hold on
subplot(3,1,2)
hold on
subplot(3,1,3)
hold on
for i = 1:simiter
    subplot(3,1,1)
    plot(simt, F1{i}, 'color', col(i,:))
    
    subplot(3,1,2)
    plot(simt, F2{i}, 'color', col(i,:))
    
    Ftendon = Fmax*k*(lt{i}-tslackl).*heaviside(lt{i}-tslackl);
    subplot(3,1,3)
    plot(simt, Ftendon, 'color', col(i,:))
    
end
subplot(3,1,1)
ylabel('F_1')
subplot(3,1,2)
ylabel('F_2')
subplot(3,1,3)
ylabel('F_{tendon}')
xlabel('Time (s)')

% Plot velocities
figure()
hold on
box on

plot(simt, v1{1}/lopt1, 'r')
plot(simt, v2{1}/lopt2, 'b')
xlabel('Time (s)')
ylabel('Velocity (Lopt/s)')


% Plot angles
figure()
subplot(2,1,1)
hold on
subplot(2,1,2)
hold on
for i = 1:simiter
    subplot(2,1,1)
    plot(simt, angle1{i}, '.-')
    subplot(2,1,2)
    plot(simt, angle2{i}, '.-')
end
ylabel('angle')