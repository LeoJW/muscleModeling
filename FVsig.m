function [FV] = FVsig(s,v)
% Sigmoid version of FV function

s1 = s(1);
s2 = s(2); 
s3 = s(3); % affects steepness of slope at 0
cmax = s(4); % "asymptote", upper limit
s4 = 1/cmax;
vmax = s(5);

v = v./vmax;
FV = s1./(s4 + s2.*exp(-s3*v));

end

