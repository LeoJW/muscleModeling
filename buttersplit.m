function [flap] = buttersplit(tdata,rawdata)
% Function for cleaning up data and splitting + overlaying waveforms
% tdata is time vector from raw data
% tnew is desired new time vector
%--Trim data
convertNaN = isnan(rawdata);
dataTrimmed = rawdata(find(convertNaN<1,1):find(convertNaN<1,1,'last'));
timeTrimmed = tdata(find(convertNaN<1,1):find(convertNaN<1,1,'last'));
%--Interpolate gaps in data
dataFilled = spline(timeTrimmed,dataTrimmed,timeTrimmed);
%--Apply LPF
[beep,boop] = butter(5,0.3);
dataFilt = filtfilt(beep,boop,dataFilled);
%--Find local minima w/ findpeaks
[~,locs] = findpeaks(-dataFilt);
%---Preallocate for loop
nwaves = length(locs)-1;
tVecNew = NaN(1,length(timeTrimmed));
for i = 1:nwaves
    % Convert time vector so all waves go from 0 to 1 (dimensionless)
    tVecNew(locs(i):locs(i+1)) = linspace(0,1,locs(i+1)-locs(i)+1);
end
% Trim NaNs out

tVecNew(isnan(tVecNew)) = [];
dataIn = dataFilt(locs(1):locs(end));
% Smooth overlaid cycles
[tOut,dataOut] = prepareCurveData(tVecNew,dataIn);
flap = fit(tOut,dataOut,'smoothingspline','SmoothingParam',0.995);
%flap = smooth(tOut,dataOut,'lowess');
end

