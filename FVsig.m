function [FV] = FVsig(s,v)
% Sigmoid version of FV function

s1 = s(1);
s2 = s(2); % 1/cmax, "asymptote", upper limit
s3 = s(3);
s4 = s(4); % affects steepness of slope at 0

v = v./vmax;
FV = s1./(s2 + s3.*exp(-s4*v));

end

