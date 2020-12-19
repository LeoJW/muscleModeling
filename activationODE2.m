function [a] = activationODE2(u,d,gam1,gam2,fsamp)
%d - delay time (assumed to be in ms)

beta1 = gam1 + gam2;
beta2 = gam1.*gam2;
alpha = (1 + beta1 + beta2);

a = zeros(size(u));

% Figure out how many samples to go back roughly 1ms
if fsamp <= 1e3 % case when each sample is larger than 1 ms
    ind = 1;
else %otherwise round to whatever is roughly 1 ms
    ind = round(1*fsamp/1e4); 
end

%a(1:d) = u(1:d); this line causes issues with first cycle of activation

for i = (d+1):length(u)
    a(i) = alpha*u(i-d) - beta1*a(i-ind) - beta2*a(i-2*ind);
end

a(a > 1) = 1;
end

