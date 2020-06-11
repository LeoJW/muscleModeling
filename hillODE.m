function [h] = hillODE(x,F)
% ODE for MTU
x = ;
F = ;
l = ;
xdot = ; % dx/dt
Fdot = ; % dF/dt
ldot = ; % dl/dt

FLactFunc = @(b,x) exp(-(((x-b(2))-1)./b(1)).^2);
b1 = 0.25; % FLact
b2 = 0; % FLact
b = [b1,b2];
FLact = FLactFunc(b,x);

p1 = 4; % FLpas
p2 = 1; % FLpas
p = [p1,p2];
FLpas = FLpasFunc(p,x)

% FV sigmoid function
Fmax = 1; % maximum force in N
cmax = 1.8; % asymptote as v approaches -inf
s = [1 1/cmax];
g = 6; % affects steepness of slope at 0
z = 0.5;

FVsig = s(1)./(s(2) + z.*exp(-g*v));

xdot = (ln(((s(1)*Fmax.*FLact.*a)/(k(l-x)-FLpas) - s(2))./z))./g;
Fdot = k(ldot - xdot);

end

%% More functions

%Passive force-length curve function
function [y] = FLpasFunc(p,x)
    y = p(1).*(x-p(2)).^2;
    y(x<p(2)) = 0;
end
