% KINEMATICS-DRIVEN EMR MODEL SCRIPT
clear vars; close all;


%% Constants




%% Define functions

FLactFunc = @(c,x) exp(-(((x-c(2))-1)./c(1)).^2);


%% Load data

kinematicsData = readtable("2019_07_02_Tae_gut_allFlights.csv");

% Assign variables for coordinates of each point
humerus_x = kinematicsData.humerus_x;
humerus_y = kinematicsData.humerus_y;
humerus_z = kinematicsData.humerus_z;
elbow_x = kinematicsData.elbow_x;
elbow_y = kinematicsData.elbow_y;
elbow_z = kinematicsData.elbow_z;
wrist_x = kinematicsData.wrist_x;
wrist_y = kinematicsData.wrist_y;
wrist_z = kinematicsData.wrist_z;
cmc_x = kinematicsData.cmc_x;
cmc_y = kinematicsData.cmc_y;
cmc_z = kinematicsData.cmc_z;

% Create n*3 matrices for each point
humerus = [humerus_x, humerus_y, humerus_z];
elbow = [elbow_x, elbow_y, elbow_z];
wrist = [wrist_x, wrist_y, wrist_z];
cmc = [cmc_x, cmc_y, cmc_z];

% Create vectors for distances between points
humerus_elbow_x = humerus_x - elbow_x;
humerus_elbow_y = humerus_y - elbow_y;
humerus_elbow_z = humerus_z - elbow_z;

humerus_elbow = [humerus_elbow_x, humerus_elbow_y, humerus_elbow_z];

elbow_wrist_x = elbow_x - wrist_x;
elbow_wrist_y = elbow_y - wrist_y;
elbow_wrist_z = elbow_z - wrist_z;

elbow_wrist = [elbow_wrist_x, elbow_wrist_y, elbow_wrist_z];

% Calculate angle between two vectors
ab = dot(humerus_elbow, elbow_wrist);
mag_a = sqrt((humerus_elbow_x).^2 + (humerus_elbow_y).^2 + (humerus_elbow_z).^2);
mag_b = sqrt((elbow_wrist_x).^2 + (elbow_wrist_y).^2 + (elbow_wrist_z).^2);

cosTheta_elbow = ab./(mag_a.*mag_b);
Theta_elbow = acos(cosTheta_elbow);

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
