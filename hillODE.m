function [dsdt] = hillODE(t,tvec,l,ldot,xF,a,C2)
% ODE for MTU

% Constants
b1 = C2(1);
b2 = C2(2);
p1 = C2(3);
p2 = C2(4);
s1 = C2(5);
s2 = C2(6);
s3 = C2(7);
s4 = C2(8);
Fmax = C2(9);

b = C2(1:2);
p = C2(3:4);
s = C2(5:8);

% FL active component
FLactFunc = @(b,x) exp(-(((x-b(2))-1)./b(1)).^2);
FLact = FLactFunc(b,x);

% FL passive component
    function [y] = FLpasFunc(p,x)
    y = p(1).*(x-p(2)).^2;
    y(x<p(2)) = 0;
    end
FLpas = FLpasFunc(p,x);

% FV sigmoid function for ref: FVsig = s1/(s2 + s3.*exp(s4*v));
% rearranged below to solve for xdot (velocity)

tvec = linspace(1,10,1e4)

x = interp1(xt,x,t);
F = interp1(Ft,F,t);
xF = [x,F];

dsdt = zeros(2,1) % [xdot, Fdot]
dsdt(1) = (ln(((s1*Fmax.*FLact.*a)/(k(l-x)-FLpas) - s2)./s3))./s4; % dx/dt
dsdt(2) = k(ldot - (ln(((s1*Fmax.*exp(-((((l-(F/k))-b2)-1)./b1).^2).*a)/(F-p1.*((l-(F/k))-p2).^2) - s2)./s3))./s4); % dF/dt

end

