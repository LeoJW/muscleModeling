%% READ KINEMATICS AND EMG FROM MONTANA STARLING WIND TUNNEL DATA
function [kine,emg] = readKineEMG()




%% Previously Known Info
%My solution to a few problems of synchronizing and loading the EMG and
%kinematics data requires the program know a few things about the data it
%can only get from a user figuring it out. All that data is stored here

%---Basic Constants
%Sampling rates
fsamp = 10e3;  %[Hz] EMG data
fksamp = 250; %[Hz] kinematics data
%Voltage Tolerance below which to count camera trigger 
trigTol = -3;
%Expected span of trigger pulse, in samples
trigSpan = 5000; 
%Max number of runs any bird has
numRuns = 13;

%---Note what last frame is for this video
lastFrame = -215; %Bird 3 run 07



%% ESTABLISH FILE PATHS
%Note initial directory
startdir = pwd;
addpath(startdir)
%Note the data directories
emgdir = strcat(startdir,'/StarlingEMG/');
kinedir = strcat(startdir,'/StarlingKinematics/');


%% LOAD DATA
%Save data into cell array for each channel (BIC, TPB, TPL, PEC, TRIG)
%Save in format XX{bird#, run#}

%---LOAD EMG DATA---%
%Preallocate (could be faster but whatever)
bic = [];
tpb = [];
tpl = [];
pec = [];
trg = [];
t = [];

%Jump to EMG data dir
cd(emgdir)
%Get info on folders within
datadirInfo = dir;
datadirInfo(1:3) = [];

%Jump to the only dir around
cd(strcat(emgdir,'/',datadirInfo(1).name))
%Get info on this dir
thisdir = dir;
%Pull names of runs from this dir
thisdirNames = {thisdir.name};
%Take only .mat files
thisdirNames(~contains(thisdirNames,'.mat')) = [];


%---LOAD THAT DATA---%
load(thisdirNames{1})

%Remove empty channel (happens when there's 6 channels)
if length(dataend) == 6
    datastart(3,:) = [];
    dataend(3,:) = [];
end
datastart = datastart(:,1);
dataend = dataend(:,1);

%Save actual traces
bic = data(datastart(1):dataend(1));
tpb = data(datastart(2):dataend(2));
tpl = data(datastart(3):dataend(3));
pec = data(datastart(4):dataend(4));
trg = data(datastart(5):dataend(5));
t = (1:length(bic))/fsamp;



%---LOAD KINEMATICS DATA---%
%Jump to kinematics directory
cd(kinedir)
%Get info on folders within
kinedirInfo = dir;
kinedirInfo(1:3) = [];
kinedirInfo([kinedirInfo.isdir]~=1) = []; %remove everything but folders

%Loop through each date folder (not really jk)
i = 1;
%Pop into date dir
cd(strcat(kinedir,kinedirInfo(i).name,'/'))
%Get info on files in this dir
datedir = dir;
ptsdirnames = {datedir.name};

%---Find csv files with points
%Keep only .csv files
ptsdirnames(~contains(ptsdirnames,'.csv')) = [];
%Keep only file with xyz points saved
ptsdirnames(~contains(ptsdirnames,'xyzpts')) = [];

%---Save XYZ
xyz = dlmread(ptsdirnames{1},',',1,0);

%Note length of kinematics data for creating time vector later
kineL = length(xyz);

%Jump back to original directory
cd(startdir)

%% Synchronize EMG and Kinematics Data

%Preallocate kinematics arrays based on lengths of data
tk = zeros(length(xyz),1);
theta = NaN(length(xyz),1);
phi = NaN(length(xyz),1);


%---Pick out camera triggers (assuming at most 2 trigs)---%
trgDiff = diff(trg);
%First trigger edge
[tv1,ti1] = min(trgDiff); %Get edge
%Remove all values around edge so not picked up again
trgDiff(ti1-trigSpan:ti1+trigSpan) = Inf; 
%Second trigger edge
[tv2,ti2] = min(trgDiff); %Get edge

%Find which of the 2 are real trigs
realTrig = [tv1,tv2] < trigTol;
%If there's more than one real trig edge & later one somehow comes first
if all(realTrig) && ti1 > ti2
    %Put the later one at position 2
    [ti1,ti2] = swap(ti1,ti2);
%If only real trig comes first, swap so real trig is in 2nd position
elseif realTrig(1) && ~realTrig(2)
    [ti1,ti2] = swap(ti1,ti2);
%Else if no real trigs (Bird5 Run01)
elseif all(~realTrig)
    ti1 = 1;
    ti2 = length(trg);
end
%Save to main matrices
trgInd = [ti1,ti2]+1; %+1 to correct diff function
%Cut all past last camera trigger
bic(trgInd(2):end) = [];
tpb(trgInd(2):end) = [];
tpl(trgInd(2):end) = [];
pec(trgInd(2):end) = [];
t(trgInd(2):end) = [];
%Zero time at trigger
t = t - t(end);



%---Create time vectors for kinematics data---%
%Define kinematics time vector using known last frame, sample rate
tk = (-kineL + (lastFrame))/fksamp :1/fksamp: (lastFrame-1)/fksamp;
tk = tk + 1/fksamp;

%---Pick out and snip kinematics data---%
%Keep only EMG data overlapping kinematics
%Find where EMG data overlaps kine data
[~,ffi] = min(abs(t-tk(1))); %First frame index
[~,lfi] = min(abs(t-tk(end))); %Last frame index
%Generate indices of EMG to keep
snipInds = ffi:lfi;
%Rewrite original data
bic = bic(snipInds);
tpb = tpb(snipInds);
tpl = tpl(snipInds);
pec = pec(snipInds);
t = t(snipInds);


%% Processing of Data

%---Convert xyz positions into angles 
%NOTE-- xyz is structured in the form: [p1_X, p1_Y, p1_Z, p2_X, p2_Y, ...] 
%       Columns are x,y,z coordinates for each tracking point, cycling
%       through all tracking points. 
%    P1 is wrist, p2 elbow, p3 shoulder, and p4 is manus
%Elbow angle
theta = xyz2theta(xyz(:,1:9));
%Manus angle
phi = xyz2theta(xyz(:,[1:6,10:12]));





%% Group Data all together for output
%NOTE-- I can make this configurable and controllable from the outside for
%       convenience so you can choose what parts of the data you get, but
%       too lazy for that right now

kine = [tk.',theta,phi];

emg = [t,bic,tpb,tpl,pec];

end