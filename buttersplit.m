function meanClean = buttersplit(tdata,rawdata,tnew)
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
twave = zeros(1,nwaves);
for i = 1:nwaves+1
    % Define each wave/cycle duration
    wavedur(i) = locs(i+1)-locs(i);
    % Make new time vector so all waves go from 0 to 1 (dimensionless)
    twave(i) = linspace(0,1,wavedur(i));
    % Prep wave var
    wave = cell(length(???),nwaves);
    % Loop through different cycles
    for j = 1:nwaves
        % Declare vars
        wave{j} = zeros(1,twave(i));
        % Define each wave
        wave{j} = dataFilt(locs(i):locs(i+1));
    end
end
t = [twave(1:end)];
%--Take avg/overlay of all cycles
meanClean = mean(wave); % use another method
% Overlay, then use Lowess model or equiv
end

