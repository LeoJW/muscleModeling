function [t,flap] = buttersplit(tdata,rawdata,tnew)
% Function for cleaning up data and splitting + overlaying waveforms
% tdata is time vector from raw data
% tnew is desired new time vector
%--Trim data
convertNaN = isnan(rawdata);
dataTrimmed = rawdata(find(convertNaN<1,1):find(convertNaN<1,1,'last'));
timeTrimmed = tdata(find(convertNaN<1,1):find(convertNaN<1,1,'last'));
%--Interpolate gaps in data
dataFilled = spline(timeTrimmed,dataTrimmed,tnew);
%--Apply LPF
[beep,boop] = butter(5,0.3);
dataFilt = filtfilt(beep,boop,dataFilled);
%--Find local minima w/ findpeaks
[~,locs] = findpeaks(-dataFilt);
%--Use locations of minima to define beginning/end of each cycle
nwaves = length(locs)-1;
wavedur = zeros(1,nwaves);
%twave = zeros(1,nwaves);
tVecNew = zeros(1,length(timeTrimmed));
for i = 1:nwaves
    % Make new time vector so all waves go from 0 to 1 (dimensionless)
    tVecNew(locs(i):locs(i+1)) = linspace(0,1,locs(i+1)-locs(1));
end
% Smooth overlaid cycles
flap = fit([tVecNew dataFilt],'lowess');
end

