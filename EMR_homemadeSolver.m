%% EMR MTU Model: Homemade Solver Version
clear vars;

% FL curves are normalized force plotted against L/Lopt
% velocity in (L/Lopt)/sec normalized to v/vmax
% FV curve is normalized force plotted against v/vmax


%% Declare comstants, setup

%---Primary controls
simiter = 8; %# of spring constants to compare
h = 1e-3; %simulation step size
velBruteForceSize = 1e4; %# of points to solve for velocity at

k = linspace(0.01,0.2,simiter); % Spring constant

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
startdur = ceil(tstart*lcycle); % start of activation in cycle
enddur = round(startdur + duration*lcycle); % duration of cycle activated in 1/1e4 s

b1 = 0.25; % FLact
b2 = 0; % FLact
b = [b1,b2];
p1 = 4; % FLpas
p2 = 1; % FLpas
p = [p1,p2];

Fmax = 1; % maximum force in N
cmax = 1.8; % asymptote as v approaches -inf
vmax = 1; % maximum velocity within range Wakeling (2012), Josephson (1993)
c1 = 0.29; % from Biewener et al. (2014)
c2 = 1; % overall curvature of FV
fvc = [c1,c2,cmax,vmax];

delay = 50; % activation delay, in ms -> need to rescale in a
gam1 = -0.993; % activation constant
gam2 = -0.993; % activation constant

% FV curve sigmoid version
% cmax same as above
s1 = 1.8;
s4 = 1;
s2 = s1-s4; % Enforces FV(0)=1
s3 = 6; % affects steepness of slope at 0
% s = [s1,s2,s3,cmax,vmax];

% Neural excitation, vector of zeros except one chunk which is 1s
ucycle = zeros(1,lcycle);
ucycle(startdur:enddur) = 1;
u = repmat(ucycle,1,ncycles);

% Activation function
d = (delay*1e-3)*(niter/totaltime); % delay, scaled
a = activationODE2(u,d,gam1,gam2);


%Singularity adjustments
Ftol = 0.1; %Tolerance for F to avoid singularities
atol = 0.01; %Tolerance for a to avoid singularities
a = (1-atol).*a+atol;


%---Anonymous functions
%FL active component function
FLactFunc = @(b,x) exp(-(((x-b(2))-1)./b(1)).^2);
%FL passive component function
FLpasFunc = @(p,x) heaviside(x-p(2)).*p(1).*(x-p(2)).^2;


%% Run Simulation

% Vector input for hillODE constants
C = [b1,b2,p1,p2,s1,s2,s3,s4,Fmax,k];

%---MTU Overall length/velocity parameters
wr = 2*pi*w; % radians per second
Lamplitude = 0.2; % amplitude of l
l = Lamplitude.*sin(wr.*t) + 2; % MTU length, l/Lopt
ldot = Lamplitude.*wr.*cos(wr*t); % MTU velocity, ldot/vmax


%Prepare Figure for loop
close all
figure(1)
hold on
box on
grid on
col = copper(simiter);

%Create simulation time vector
simt = 0:h:totaltime;
%Prepare variables for loop
err = cell(size(k));
Ferr = cell(size(k));
v = cell(size(k));
x = cell(size(k));
F = cell(size(k));
vsweep = linspace(-1,1,velBruteForceSize);
wrk = cell(size(k));
pwr = cell(size(k));
%Calculate FV function at all velocities
FVactVal = FVsig([s1,s2,s3,s4,vmax],vsweep);

%Loop through different spring constants
for i = 1:simiter
    
    %Declare vectors for this run's simulation
    x{i} = [1,zeros(1,length(simt)-1)]; %initial condition for muscle length
    v{i} = zeros(1,length(simt));
    err{i} = zeros(1,length(simt));
    F{i} = zeros(1,length(simt));
    wrk{i} = zeros(1,length(simt));
    pwr{i} = zeros(1,length(simt));
    
    %Loop thru each point in time for simulation
    for j = 1:length(simt)
        %Find new x from previous velocity
        if j~=1
            x{i}(j) = x{i}(j-1) + v{i}(j-1)*h;
        end
        %interpolate l,a at this time
        tl = interp1(t,l,simt(j));
        ta = interp1(t,a,simt(j));
        %Solve individual components of hill model
        FLactVal = (1-Ftol).*FLactFunc([b1,b2],x{i}(j)) + Ftol;
        %use x(j) to solve for muscle velocity
        eval = k(i)*(tl-x{i}(j)) - (FLactVal.*FVactVal.*ta + FLpasFunc([p1,p2],x{i}(j)));
        %Find root of function where velocity is valid
        [errval,vind] = min(abs(eval));
        err{i}(j) = eval(vind);
        v{i}(j) = vsweep(vind);
        F{i}(j) = hillv2(x{i}(j),v{i}(j),ta,B);
        
        % work, area under curve w/ neg vs pos velocity
        wrk{i}(j) = trapz(x{i}(j),F{i}(j).*sign(v{i}(j)));
        % instantaneous power
        pwr{i}(j) = F{i}(j).*v{i}(j);
        
    end

    %Plot output
    plot(simt, v{i},'color',col(i,:))
    drawnow
    
end

%Aesthetics
xlabel('Time (s)')
ylabel('Muscle Velocity')
%Aesthetics for colorbar
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
