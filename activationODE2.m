function [a] = activationODE2(u,d,gam1,gam2)
%d - delay time (assumed to be in ms)

beta1 = gam1 + gam2;
beta2 = gam1.*gam2;
alpha = (1 + beta1 + beta2);

a = zeros(size(u));

a(1:d) = u(1:d);

for i = (d+1):length(u)
    a(i) = alpha*u(i-d) - beta1*a(i-1) - beta2*a(i-2);
end

% a(a > 1) = 1;
end

