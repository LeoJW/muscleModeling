%% Practice with simple mass-spring system

% Constants
k = 0.5; % spring constant
M = 5; % mass
niter = 1000 % number of repeats in for loop
xi = 2; % initial position
vi = 5; % initial velocity
ai = 8; % initial acceleration

% Variables
t = linspace(1,100,niter)
x(t) = x(t-1) + v(t-1).*t;
a(t) = (k./M).*x(t-1);
v(t) = x(t-1) + a(t-1).*t;

% Make loop for spring motion

figure()

for i = 1:niter
    y1 = x(t);
    plot(t,y1,'color',col(i,:));
end
