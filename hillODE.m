function [xdot,Fdot] = hillODE(t,l,ldot,x,F,a)
% ODE for MTU

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

xt = linspace(1,10,1e4);
Ft = linspace(1,10,1e4);

x = interp1(xt,x,t);
F = interp1(Ft,F,t);

xdot = (ln(((s1*Fmax.*FLact.*a)/(k(l-x)-FLpas) - s2)./s3))./s4; % dx/dt
Fdot = k(ldot - (ln(((s1*Fmax.*exp(-((((l-(F/k))-b2)-1)./b1).^2).*a)/(F-p1.*((l-(F/k))-p2).^2) - s2)./s3))./s4); % dF/dt

end

