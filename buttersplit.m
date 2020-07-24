function meanClean = buttersplit(tdata,rawdata,tnew)
% Function for cleaning up data and splitting + overlaying waveforms
%--Trim data
convertNaN = isnan(rawdata);
dataTrimmed = rawdata(find(convertNaN>0,1):find(convertNaN>0,1,'last'));
%--Interpolate gaps in data
dataFilled = spline(tdata,dataTrimmed,tnew);
%--Apply LPF
[beep,boop] = butter(5,0.5);
dataFilt = filtfilt(beep,boop,dataFilled);
%--Find local minima w/ findpeaks
[~,locs] = findpeaks(-dataFilt);
%--Use locations of minima to define beginning/end of each cycle
nwaves = length(locs)+1;
for i = 1:nwaves
    % Define each wave/cycle duration
    wavedur(i) = length(locs(i):locs(i+1));
    % Prep wave var for next loop
    wave = cell(wavedur(i),nwaves);
    % Loop through different cycles
    for j = 1:length(wavedur(i))
        % Declare vars
        wave{j} = zeros(1,wavedur(i));
        % Define each wave
        wave{j} = rawdata(locs(i):locs(i+1));
    end
end
%--Take avg of all split cycles
meanClean = mean(wave);
end

