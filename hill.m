% Hill function
function [H] = hill(x,v,a,C)

b1 = C(1);
b2 = C(2);
p1 = C(3);
p2 = C(4);
c1 = C(5);
c2 = C(6);
cmax = C(7);
vmax = C(8);
u = C(9);
d = C(10);
gam1 = C(11);
gam2 = C(12);
Fmax = C(13);

FLactFunc = @(b,x) exp(-(((x-b(2))-1)./b(1)).^2); % FL active component

    function [y] = FLpasFunc(p,x) % FL passive component
    y = p(1).*(x-p(2)).^2;
    y(x<p(2)) = 0;
    end

F = FV4param(c,v); % Force velocity curve
c1 = c(1); %Curve shape during contracting phase
c2 = c(2); %Curvature
c3 = (c1.^2 + c1)/c1.^2; %Slope at zero 
cmax = c(3); %Asymptote as v approaches -inf
vmax = c(4); %Max velocity; assumes velocity input is actual length/s
kl = (c2*(1-cmax))/(c3*(1+c2));

v = v./vmax;
F = zeros(size(v));

for i = 1:length(v)
    if v(i) <= 0
        F(i) = (1 + v(i))/(1 - v(i)/c1);
    else
        F(i) = (kl - cmax*v(i))/(kl-v(i));
    end
end

a = activationODE2(u,d,gam1,gam2); % Activation function
%d - delay time (assumed to be in ms)

beta1 = gam1 + gam2;
beta2 = gam1.*gam2;
alpha = (1 + beta1 + beta2);

a = zeros(size(u));

a(1:d) = u(1:d);

for i = (d+1):length(u)
    a(i) = alpha*u(i-d) - beta1*a(i-1) - beta2*a(i-2);
end

% a(a > 1) = 1;

H = Fmax.*((FLactFunc.*F.*a) + y);
end