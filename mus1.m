function [wrk] = mus1(control,w,C,conv,sing)
% Function version of EMRscript.m


%% Controls and constants

%--Sim controls
simiter = control(1); % number of activation phases to compare
h = control(2); % step size
velBruteSize = control(3); % number of points to solve for v
stimPhase = control(4); % version of tstart that varies
ncycles = control(5); % number of cycles
duration = control(6); % activ dur

%--Hill constants
b1 = C(1); % FLact
b2 = C(2); % FLact
p1 = C(3); % FLpas
p2 = C(4); % FLpas

cmax = C(5); % asymptote as v approaches -inf
vmax = C(6); % maximum velocity in Lopt/s
c1 = C(7); % from Biewener et al. (2014)
c2 = C(8); % overall curvature of FV

b = C(1:2);
p = C(3:4);
fvc = C(5:8);

m1 = C(9); % scaling factor
m2 = C(10); % horizontal translation
m3 = C(11); % slope
m = C(9:11); % FV curve, smooth ramp portion
delay = C(12); % activation delay, in ms -> rescaled in a
gam1 = C(13); % activation constant
gam2 = C(14); % activation constant

%--Conversion constants
mRL = conv(1); % muscle resting length, mm
mtuRL = conv(2); % mtu resting length, mm
Lopt = mRL + 0.05*mRL; % mm
lamplitude = conv(3);
EMRArea = conv(4); % (mm^2) dry density in g/mm^3, mass in g
Fmax = 300e3*1e-6*EMRArea; % max force in N (convert from 300kPa to N/mm^2, multiply by EMR area)

tslackl = conv(5); % tendon slack length, mm
kActual = conv(6); % spring constant, N/mm
k = kActual*(Lopt/Fmax); % dimensionless

%--Singularity adjustments
Ftol = sing(1); % tolerance for F to avoid singularities
atol = sing(2); % tolerance for a to avoid singularities


%% Simulation constants setup

%--Create simulation time vector
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


%% Run simulation

%--MTU length parameters
wr = 2*pi*w; % convert to radians
l = (lamplitude.*sin(wr.*simt) + mtuRL)/Lopt; % MTU length in Lopt/s

%--Split cycles
cycNum = repelem(1:ncycles,lcycle);

% Prepare variables for loop
err = cell(1,simiter);
F = cell(1,simiter);
v = cell(1,simiter);
x = cell(1,simiter);
vsweep = linspace(-vmax,vmax,velBruteSize);
wrk = cell(1,simiter);
pwr = cell(1,simiter);
% Calculate FV function at all velocities
FVactVal = FV4param(fvc,vsweep) + FVactHinge(m,vsweep);


tic
%---Loop through different stimulation phases
parfor i = 1:simiter
    
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
    x0 = x0sweep(ind); % dimensionless, L/Lopt
    
    % Declare vectors for simulation run
    x{i} = [x0,zeros(1,length(simt)-1)]; % muscle length initial condition
    v{i} = [v0,zeros(1,length(simt)-1)]; % velocity initial condition
    err{i} = zeros(1,length(simt)); % error
    F{i} = zeros(1,length(simt)); % force
    wrk{i} = zeros(1,length(simt));
    pwr{i} = zeros(1,length(simt));
    
    % Loop through each time point
    for j = 2:length(simt)
        % Find new x from previous v
        x{i}(j) = x{i}(j-1) + v{i}(j-1)*h;
        
        % Solve individual components of Hill model
        FLactVal = (1-Ftol).*FLactFunc([b1,b2],x{i}(j)) + Ftol;
        % Use x(j) to solve for muscle v
        eval = k*(l(j)-x{i}(j)-tslackl/Lopt).*heaviside(l(j)-x{i}(j)-tslackl/Lopt) - ...
            (FLactVal.*FVactVal.*a{i}(j) + FLpasFunc([p1,p2],x{i}(j)));
        % Find root of function where velocity is valid
        [errval,vind] = min(abs(eval));
        err{i}(j) = eval(vind);
        v{i}(j) = vsweep(vind);
        F{i}(j) = hill(x{i}(j), v{i}(j), a{i}(j), C);
        
        % Re-find x from currently calculated v (backwards euler)
        x{i}(j) = x{i}(j-1) + h*v{i}(j);
        
        % instantaneous power
        pwr{i}(j) = F{i}(j).*v{i}(j);
        
    end
    
    % Convert values to real units
    x{i} = x{i}*Lopt; % converts length to mm
    F{i} = F{i}*Fmax; % converts force to N
    
    % work, area under curve w/ neg vs pos velocity
    wrk{i} = -trapz(x{i}(cycNum>(ncycles-1)), F{i}(cycNum>(ncycles-1))); % units of N*mm, or mJ
    
end

toc

end

