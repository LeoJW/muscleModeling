function [dxdt] = hillODE(t,tvec,l,ldot,state,a,C)
% ODE for MTU

% State vector
x = state(1);
F = state(2);

% Constants
b1 = C(1);
b2 = C(2);
p1 = C(3);
p2 = C(4);
s1 = C(5);
s2 = C(6);
s3 = C(7);
s4 = C(8);
Fmax = C(9);
k = C(10);

b = C(1:2);
p = C(3:4);

% FL active component
FLactFunc = @(b,x) exp(-(((x-b(2))-1)./b(1)).^2);
FLact = FLactFunc(b,x);

% FL passive component
    function [y] = FLpasFunc(p,x)
    y = p(1).*(x-p(2)).^2;
    y(x<p(2)) = 0;
    end
FLpas = FLpasFunc(p,x);

% could pass in FLpas and FLact as arguments to hillODE instead
% FV sigmoid function for ref: FVsig = s1/(s2 + s3.*exp(s4*v));
% rearranged below to solve for xdot (velocity)

% l and ldot at specific time points
tl = interp1(tvec,l,t);
tldot = interp1(tvec,ldot,t);
ta = interp1(tvec,a,t);

dxdt = zeros(2,1); % [xdot, Fdot]
if a > 0
    dxdt(1) = (-1/s3).*(log(((s1/s2)*(Fmax.*FLact.*ta)/(k*(tl-x)-FLpas) - s4/s2))); % dx/dt
    dxdt(2) = k*(tldot - ((-1/s3).*log(((s1/s2)*(Fmax.*exp(-((((tl-F/k)-b2)-1)./b1).^2).*ta)/(F-p1.*((tl-F/k)-p2).^2) - s4/s2)))); % dF/dt
else
    dxdt(1) = tldot./((H(x-p2).*(2p1.*(x-p2).^2))/k + 1))
    dxdt(2) = k*(tldot - tldot./((H((tl-F/k)-p2).*(2*p1*((tl-F/k)-p2).^2))/k + 1))

end

