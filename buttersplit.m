function [outputArg1,outputArg2] = buttersplit(tdata,rawdata,tnew)
% Function for applying LPF to data and splitting waveforms
% Trim data
convertNaN = isnan(rawdata);
dataTrimmed = rawdata(find(convertNaN>0,1):find(convertNaN,'last'));
% Interpolate gaps in data
dataFull = spline(tdata,dataTrimmed,tnew);
% Apply LPF
[beep,boop] = butter(5,0.5);
dataFilt = filtfilt(beep,boop,dataFull);
% Find local minima w/ findpeaks
[mins,locs] = findpeaks(-dataFilt);
% Use locations of minima to define beginning/end of each cycle
wave = zeros(1,length(locs));
for i = 1:length(locs)
    wave(i) = thetaraw(locs(i):locs(i+1));
end
end

