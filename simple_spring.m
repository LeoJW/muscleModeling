%% Practice with simple mass-spring system

% Constants
k = 0.5; % spring constant
M = 5; % mass
t(1) = 0; % start time
x(1) = 2; % initial position
v(1) = 0; % initial velocity
dt = 1 % time step

% Make loop for spring motion

niter = 1000
t = zeros(1,niter+1);
x = zeros(1,niter+1);
v = zeros(1,niter+1);
a = zeros(1,niter+1);

for i = 1:niter
    t(i+1) = t(i) + dt;
    a(i+1) = (k./M).*x(i);
    v(i+1) = v(i) + a(i).*dt;
    x(i+1) = x(i) + v(i).*dt;
end

plot(t,x)
xlabel("Time")
ylabel("Distance")
