%% Practice with simple mass-spring system
clear vars; 


% Constants
k = 0.5; % spring constant
M = 5; % mass
t(1) = 0; % start time
x(1) = 2; % initial position
v(1) = 0; % initial velocity
dt = 0.01; % time step

% Make loop for spring motion

niter = 1000;
t = zeros(1,niter+1);
x = zeros(1,niter+1);
v = zeros(1,niter+1);
a = zeros(1,niter+1);

for i = 1:niter
    t(i+1) = t(i) + dt;
    a(i+1) = -(k./M).*x(i);
    v(i+1) = v(i) + a(i).*dt;
    x(i+1) = x(i) + v(i).*dt;
end

figure(1)
plot(t,x)
xlabel("Time")
ylabel("Distance")


%% Leo's Edited Version

%Constants
k = 5; % spring constant
M = 5; % mass
simTime = 10; %seconds
dt = 0.1; % time step should be small: the smaller, the more accurate
t = 0:dt:simTime; %<-Many ways to define time (e.g. as a linspace, or like this),
%               but code will run faster when not saving/defining things in a loop
niter = length(t);


%Initial Conditions
x0 = [2,0]; %[position,velocity]


%Preallocate for loop
x = [x0(1), zeros(1,niter-1)];
v = [x0(2), zeros(1,niter-1)];
a = [(k/M)*x0(1), zeros(1,niter-1)];

for i = 2:niter
    a(i) = (-k/M)*x(i-1);
    v(i) = v(i-1) + a(i-1)*dt;
    x(i) = x(i-1) + v(i-1)*dt;
end

figure(2)
plot(t,x)
xlabel("Time")
ylabel("Distance")

%% Damped spring

%Constants
k = 5; % spring constant
b = 2; % damping constant
M = 5; % mass
simTime = 10; %seconds
dt = 0.1; % time step should be small: the smaller, the more accurate
t = 0:dt:simTime; %<-Many ways to define time (e.g. as a linspace, or like this),
%               but code will run faster when not saving/defining things in a loop
niter = length(t);


%Initial Conditions
x0 = [2,0]; %[position,velocity]


%Preallocate for loop
x = [x0(1), zeros(1,niter-1)];
v = [x0(2), zeros(1,niter-1)];
a = [((k/M)*x0(1) + (b/M)*x0(2)), zeros(1,niter-1)];

for i = 2:niter
    a(i) = (-k/M)*x(i-1) + (-b/M)*v(i-1);
    v(i) = v(i-1) + a(i-1)*dt;
    x(i) = x(i-1) + v(i-1)*dt;
end

figure(3)
plot(t,x)
xlabel("Time")
ylabel("Distance")