function [F] = FV4param(c,v)

c1 = c(1); %Curve shape during contracting phase
c2 = c(2); %Curvature
c3 = (c1.^2 + c1)/c1.^2; %Slope at zero 
cmax = c(3); %Asymptote as v approaches -inf
vmax = c(4); %Max velocity; assumes velocity input is actual length/s
kl = (c2*(1-cmax))/(c3*(1+c2));

v = v./vmax;
F = zeros(size(v));

for i = 1:length(v)
    if v(i) <= 0
        F(i) = (1 + v(i))/(1 - v(i)/c1);
    elseif v(i) >= 0.6
        F(i) = 0.15*v(i) + 1.604;
    else
        F(i) = (kl - cmax*v(i))/(kl-v(i)); 
    end
end

end




