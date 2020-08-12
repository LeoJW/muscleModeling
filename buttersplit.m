function [tOut,dataOut,cycDur] = buttersplit(tdata,rawdata,butterOrder,butterFreq)
% Function for cleaning up data and splitting + overlaying waveforms
% tdata is time vector from raw data
% tnew is desired new time vector
%--Trim data
convertNaN = isnan(rawdata);
dataTrimmed = rawdata(find(convertNaN<1,1):find(convertNaN<1,1,'last'));
timeTrimmed = tdata(find(convertNaN<1,1):find(convertNaN<1,1,'last'));
%--Interpolate gaps in data
nonans = isnan(dataTrimmed);
dataFilled = spline(timeTrimmed(~nonans),dataTrimmed(~nonans),timeTrimmed);
%--Apply LPF
[beep,boop] = butter(butterOrder,butterFreq);
dataFilt = filtfilt(beep,boop,dataFilled);
%--Find local minima w/ findpeaks
[~,locs] = findpeaks(-dataFilt);
% Find mean cycle duration
cycDur = mean(diff(locs));
%---Preallocate for loop
nwaves = length(locs)-1;
tVecNew = NaN(length(timeTrimmed),1);
for i = 1:nwaves
    % Convert time vector so all waves go from 0 to 1 (dimensionless)
    tVecNew(locs(i):locs(i+1)) = linspace(0,1,locs(i+1)-locs(i)+1);
end
% Trim NaNs out
tVecNew(isnan(tVecNew)) = [];
dataIn = dataFilt(locs(1):locs(end));
% Smooth overlaid cycles
[tOut,dataToSmooth] = prepareCurveData(tVecNew,dataIn);
dataSmooth = fit(tOut,dataToSmooth,'smoothingspline','SmoothingParam',0.995);
dataOut = dataSmooth(tOut);
end

