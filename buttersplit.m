function [dataSmooth,cycDur,w] = buttersplit(tdata,rawdata,butterOrder,butterFreq)
% Function for cleaning up data and splitting + overlaying waveforms
% tdata = time vector from raw data
% rawdata = raw y values
% butterOrder and butterFreq = parameters for Butterworth LPF
% nSamples = sample rate
% tOut = normalized time vector based on avg cycle duration
% dataOut = smoothed fit of 'average' overlaid cycles
% cycDur = the mean cycle duration
% w = cycle frequency from data in Hz

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
    % Convert time vector so all waves go from 0 to avg cycle length (s)
    tVecNew(locs(i):locs(i+1)) = linspace(0,1,locs(i+1)-locs(i)+1);
end
% Trim NaNs out
tVecNew(isnan(tVecNew)) = [];
dataIn = dataFilt(locs(1):locs(end));
% Smooth overlaid cycles
[tSmooth,dataToSmooth] = prepareCurveData(tVecNew,dataIn);
dataSmooth = fit(tSmooth,dataToSmooth,'smoothingspline','SmoothingParam',0.995);
w = round(length(tdata)/(cycDur*max(tdata)));

end

