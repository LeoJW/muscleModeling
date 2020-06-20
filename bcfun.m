function [res] = bcfun(ya,yb)
% Boundary conditions for bvp4c

res = [ya(1)-1
        yb(1)]

end

