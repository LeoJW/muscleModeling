function meanClean = buttersplit(tdata,rawdata,tnew)
% Function for cleaning up data and splitting + overlaying waveforms
%--Trim data
convertNaN = isnan(rawdata);
dataTrimmed = rawdata(find(convertNaN>0,1):find(convertNaN>0,1,'last'));
%--Interpolate gaps in data
dataFull = spline(tdata,dataTrimmed,tnew);
%--Apply LPF
[beep,boop] = butter(5,0.5);
dataFilt = filtfilt(beep,boop,dataFull);
%--Find local minima w/ findpeaks
[~,locs] = findpeaks(-dataFilt);
%--Use locations of minima to define beginning/end of each cycle
wave = zeros(1,length(locs));
for i = 1:length(locs)
    wave(i) = rawdata(locs(i):locs(i+1));
end
%--Take avg of all split cycles
meanClean = mean(wave);
end

