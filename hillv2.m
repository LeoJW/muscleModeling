function [H2] = hillv2(inputArg1,inputArg2)
% Hill function version 2 w/ sigmoid FV curve


% Variables in C, descriptions listed in functions
b1 = C(1);
b2 = C(2);
p1 = C(3);
p2 = C(4);
s1 = C(5);
s2 = C(6);
s3 = C(7);
s4 = C(8);
Fmax = C(9);

b = C(1:2);
p = C(3:4);
fvs = C(5:8);

FLactFunc = @(b,x) exp(-(((x-b(2))-1)./b(1)).^2); % FL active component

    function [y] = FLpasFunc(p,x) % FL passive component
    y = p(1).*(x-p(2)).^2;
    y(x<p(2)) = 0;
    end

H = Fmax.*((FLactFunc(b,x).*FVsig(fvs,v).*a) + FLpasFunc(p,x));

end

