function [H2] = hillv2(x,v,a,C2)
% Hill function version 2 w/ sigmoid FV curve


% Variables in C, descriptions listed in functions
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
fvs = C(5:8);

FLactFunc = @(b,x) exp(-(((x-b(2))-1)./b(1)).^2); % FL active component

    function [y] = FLpasFunc(p,x) % FL passive component
    y = p(1).*(x-p(2)).^2;
    y(x<p(2)) = 0;
    end

H2 = Fmax.*((FLactFunc(b,x).*FVsig(fvs,v).*a) + FLpasFunc(p,x));

end

