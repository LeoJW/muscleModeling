% KINEMATICS-DRIVEN EMR MODEL SCRIPT
clear vars; close all;


%% Constants




%% Define functions

FLactFunc = @(c,x) exp(-(((x-c(2))-1)./c(1)).^2);


%% Load data

kineData = readtable("2019_07_02_Tae_gut_allFlights.csv");

% Create n*3 matrices for each point
humerus = [kineData.humerus_x, kineData.humerus_y, kineData.humerus_z];
elbow = [kineData.elbow_x, kineData.elbow_y, kineData.elbow_z];
wrist = [kineData.wrist_x, kineData.wrist_y, kineData.wrist_z];
cmc = [kineData.cmc_x, kineData.cmc_y, kineData.cmc_z];

% Create vectors for distances between points
humerus_elbow = humerus-elbow;
elbow_wrist = elbow-wrist;
wrist_cmc = wrist-cmc;

% Calculate angle between two vectors
ab = dot(humerus_elbow, elbow_wrist, 2);
mag_a = vecnorm(humerus_elbow,2,2);
mag_b = vecnorm(elbow_wrist,2,2);

theta_elbow = acos(ab./(mag_a.*mag_b))*180/pi;

bc = dot(elbow_wrist, wrist_cmc, 2);
mag_c = vecnorm(wrist_cmc,2,2);

theta_wrist = acos(bc./(mag_b.*mag_c))*180/pi;

% Testing calculation of angle from two simple vectors
p = [1,2,3];
q = [4,5,6];

pq = dot(p,q);
mag_p = sqrt(1.^2 + 2.^2 + 3.^2);
mag_q = sqrt(4.^2 + 5.^2 + 6.^2);

cosTheta = pq./(mag_p.*mag_q);
Theta = acos(cosTheta); % Theta in radians
% It works here so there must be something funny going on in my
% humerus_elbow and elbow_wrist vectors


%% Define Hill model and fit constants




%% Run Hill model on data




%% Plot




%% More functions

%Passive force-length curve function
function [y] = FLpasFunc(c,x)
    y = c(1).*(x-c(2)).^2;
    y(x<c(2)) = 0;
end
