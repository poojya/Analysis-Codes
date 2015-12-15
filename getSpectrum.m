function [S1,mlogS1,dS1,t2,f2,average_S1] = getSpectrum(plotData,timeVals,specType,mtmParams,movingWin,takeLogTrial,BLMin,BLMax)


if ~exist ('takeLogTrial','var')
    takeLogTrial = 0;
end

if ~exist ('BLMin','var')
    BLMin = -0.5;
end

if ~exist ('BLMax','var')
    BLMax = 0;
end

if ~exist('mtmParams','var')
    mtmParams.Fs = 2000;
    mtmParams.tapers=[2 3]; % [1 1] case is simply stft with dpss window
    mtmParams.trialave=1;
    mtmParams.err=0;
    mtmParams.pad=-1;
end

if ~exist('movingWin','var')
    movingWin = [0.5 0.01];
end

if ~exist('specType','var')
    specType=3;
end

if takeLogTrial
    
    if specType==3 % line spectrum i.e. normal stft
        
        mtmParams.trialave = 0; % don't average spectrum across trials
        [S1,f2]=mtspectrumc(plotData',mtmParams);
        
        logS1 = conv2Log(S1); % log power for each trial
        mlogS1 = mean(logS1,2); % mean of log power across trials
        
        dS1 = [];
        t2 = [];
        
    else
        
        mtmParams.trialave = 1; % don't average spectrum across trials
        [S1,t2,f2]=mtspecgramc(plotData',movingWin,mtmParams);

        logS1 = conv2Log(S1); % log power for each trial
        trialmlogS1 = mean(logS1,3); % mean of log power across trials

        t2 = t2 + timeVals(1); % shift the t values to the actual time
        tBL = (t2>=BLMin) & (t2<=BLMax); % baseline time indices
        trialmlogS1BL = trialmlogS1(tBL,:); % baseline trial mean log power

        mlogS1BL = mean(trialmlogS1BL,1); % mean baseline log power

        % difference spectrum calculation
        dS1 = 10*(trialmlogS1 - repmat(mlogS1BL,size(trialmlogS1,1),1));
        
        
        % mean log power across trials
        mlogS1 = trialmlogS1;
    end
    
else
    
    if specType==3 % line spectrum i.e. normal stft
        
        mtmParams.trialave = 0; %  dont average spectrum across trials
        [S1,f2]=mtspectrumc(plotData,mtmParams);
        average_S1=mean(S1,2);
        mlogS1 = conv2Log(S1); % log power for each trial
        
        dS1 = [];
        t2= [];
        
    else
        
        mtmParams.trialave = 1; % average spectrum across trials
        [S1,t2,f2]=mtspecgramc(plotData,movingWin,mtmParams);
        t2 = t2 + timeVals(1); % shift the t values to the actual time

        tBL = (t2>=BLMin) & (t2<=BLMax); % baseline time indices
        S1BL = S1(tBL,:); % part of spectrum corresponding to the baseline period
        mlogS1BL = mean(conv2Log(S1BL),1); % mean log power across these time points at every frequency

        % difference spectrum calculation
        dS1 = 10*(conv2Log(S1) - repmat(mlogS1BL,size(S1,1),1)); % in dB
        
        % mean log power across trials
        mlogS1 = conv2Log(S1);
    end
    
end

end

