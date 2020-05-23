% Hill function
function [H] = hill(x,v,a,C)

% Variables in C, descriptions listed in functions
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

b = C(1:2);
p = C(3:4);
fvc = C(5:8);

FLactFunc = @(b,x) exp(-(((x-b(2))-1)./b(1)).^2); % FL active component

    function [y] = FLpasFunc(p,x) % FL passive component
    y = p(1).*(x-p(2)).^2;
    y(x<p(2)) = 0;
    end

F = FV4param(fvc,v); % Force velocity curve

H = Fmax.*((FLactFunc(b,x).*F(fvc,v).*a) + y(p,x))
end
