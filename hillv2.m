function [H2] = hillv2(x,v,a,B)
% Hill function version 2 w/ sigmoid FV curve

% Variables in C2, descriptions listed FL,FV functions
b1 = B(1);
b2 = B(2);
p1 = B(3);
p2 = B(4);
s1 = B(5);
s2 = B(6);
s3 = B(7);
s4 = B(8);
vmax = B(9)
Fmax = B(10);

b = B(1:2);
p = B(3:4);
s = B(5:9);

FLactFunc = @(b,x) exp(-(((x-b(2))-1)./b(1)).^2); % FL active component

    function [y] = FLpasFunc(p,x) % FL passive component
    y = p(1).*(x-p(2)).^2;
    y(x<p(2)) = 0;
    end

H2 = Fmax.*((FLactFunc(b,x).*FVsig(s,v).*a) + FLpasFunc(p,x));

end

