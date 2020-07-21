function [theta] = xyz2theta(XYZ,useDegrees)
%xyz2theta converts array of 3 XYZ points to single vector of angle between
%points. 
%Assumes XYZ is of form p1_X, p1_Y, p1_Z, p2_X, ... 
%Returns angle in degrees as default

if nargin < 2
    useDegrees = 'deg';
end

%Define vectors BA and BC
BA = XYZ(:,4:6) - XYZ(:,1:3);
BC = XYZ(:,4:6) - XYZ(:,7:9);


%Calculate theta from AB and BC
theta = acos( dot(BA,BC,2)./(vecnorm(BA,2,2).*vecnorm(BC,2,2)) );

if contains(useDegrees,'deg')
    theta = theta*180/pi;
end

end

