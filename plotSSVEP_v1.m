
function plotSSVEP_v1(subjectName,expDate,protocolName,folderSourceString,gridType)

%% Default MTM params


if ~exist('folderSourceString','var');   folderSourceString ='E:';      end
if ~exist('gridType','var');             gridType='EEG';                end
if ~exist('dataLog','var')
    dataLog=getdataLog(subjectName,expDate,protocolName,folderSourceString,gridType);
end

%% setting these flags removes the bad trials specific to each electrode
if ~exist('useSingleElectrodeBadTrialList','var'); useSingleElectrodeBadTrialList=1; end
if ~exist('useBaselineSpectrumFlag','var'); useBaselineSpectrumFlag=1; end
%% loads these parameters from the data log file
mtmParams.Fs = dataLog{9,2};
mtmParams.tapers=[1 1];
mtmParams.trialave=0;
mtmParams.err=0;
mtmParams.pad=-1;
%moving window length is  0.25 to get 4 hz resolution
movingWin = [0.25 0.01];
 
folderName = fullfile(folderSourceString,'data',subjectName,gridType,expDate,protocolName);

% Get folders
folderExtract = fullfile(folderName,'extractedData');
folderSegment = fullfile(folderName,'segmentedData');
folderLFP = fullfile(folderSegment,'LFP');
stimResults = loadStimResults(folderExtract);


% load LFP Information
[analogChannelsStored,timeVals] = loadlfpInfo(folderLFP);
load(fullfile(folderExtract,'parameterCombinations.mat'));

%%% time and frequency indexes of spectrum calculation :
tmin=0.25;
tmax=0.50;
findex=unique(stimResults.temporalFreq0Hz);
cindex=unique(stimResults.contrast0PC);
no_of_contrasts=length(cindex)+1;
blmax=0.0;
blmin=-0.25;
blPeriod = (timeVals>=blmin) & (timeVals<=blmax);
tERP = (timeVals>=tmin) & (timeVals<=tmax);
colorNames={'r','b'};

%% electrodes for analysis
occipital_list_right=[31];
occipital_list_left=[29];

%% main loop to obtain the SSVEP coefficients
 %% attend right 
                %%%% calculating baseline spectrum and subtracting it from the stimulus
%%%% time spectrum :

    %%%%%%%%
for i=2:3
    if i==2  %% for the first temporal frequency
        for j=1:2
            if j==1  %% attend right 
                for c=1:no_of_contrasts
                
                 %%%%%%%%%%%%%%%%%%%%%%%%   ATTEND 16 HZ ON THE RIGHT %%%%%%%%%%%%%%%   
                 [plotData_attend16hz_right,~,goodPos]=getDataSRC(c,2,1,1,3,folderName,folderLFP,occipital_list_left,useSingleElectrodeBadTrialList);
                 %%% baseline correction :
                  plotData_attend16hz_right_baselinecorrected=baselineCorrection(plotData_attend16hz_right,blPeriod);
                 %%spectrum estimation
                 [plotData_stimOn,plotData_baseline]=getDataForSpectrum(plotData_attend16hz_right_baselinecorrected,tERP,blPeriod);
                 [S1_attend16hz_right,mlogS1,dS1,t2,f2,average_S1_attend16hz_right] = getSpectrum(plotData_stimOn',timeVals,3,mtmParams,movingWin,0,blmin,blmax);
                 [S1_attend16hz_right_baseline,mlogS1,dS1,t2,f2,average_S1_attend16hz_right_baseline] = getSpectrum(plotData_baseline',timeVals,3,mtmParams,movingWin,0,blmin,blmax);
                 [S1,S1_power,no_of_trials]=getDifferenceSpectrum(S1_attend16hz_right,S1_attend16hz_right_baseline,average_S1_attend16hz_right,average_S1_attend16hz_right_baseline,useBaselineSpectrumFlag,findex,f2);
                 clear plotData_stimOn; clear plotData_baseline;
                 S1_attend16hz_right_spectrum{c}=S1;
                 S1_attend16hz_right_power(c)=S1_power;
                 S1_attend16hz_right_trials(c)=no_of_trials;
                 clear S1; clear S1_power; clear no_of_trials;
                 %%%%%%%%%%%%%%%%%%%%%%%%   UNATTEND 16 HZ ON THE RIGHT %%%%%%%%%%%%%%% 
                  [plotData_unattend16hz_right,~,goodPos]=getDataSRC(c,2,1,2,3,folderName,folderLFP,occipital_list_left,useSingleElectrodeBadTrialList);
                 %%% baseline correction :
                  plotData_unattend16hz_right_baselinecorrected=baselineCorrection(plotData_unattend16hz_right,blPeriod);
                 %%spectrum estimation
                 [plotData_stimOn,plotData_baseline]=getDataForSpectrum(plotData_unattend16hz_right_baselinecorrected,tERP,blPeriod);
                 [S1_unattend16hz_right,mlogS1,dS1,t2,f2,average_S1_unattend16hz_right] = getSpectrum(plotData_stimOn',timeVals,3,mtmParams,movingWin,0,blmin,blmax);
                 [S1_unattend16hz_right_baseline,mlogS1,dS1,t2,f2,average_S1_unattend16hz_right_baseline] = getSpectrum(plotData_baseline',timeVals,3,mtmParams,movingWin,0,blmin,blmax);
                 [S1,S1_power,no_of_trials]=getDifferenceSpectrum(S1_unattend16hz_right,S1_unattend16hz_right_baseline,average_S1_unattend16hz_right,average_S1_unattend16hz_right_baseline,useBaselineSpectrumFlag,findex,f2);
                 clear plotData_stimOn; clear plotData_baseline;
                 S1_unattend16hz_right_spectrum{c}=S1;
                 S1_unattend16hz_right_power(c)=S1_power;
                 S1_unattend16hz_right_trials(c)=no_of_trials;
                 clear S1; clear S1_power; clear no_of_trials;
                 
                end
            else %% attend left-first temporal frequency
               for c=1:no_of_contrasts
                
                 %%%%%%%%%%%%%%%%%%%%%%%%   ATTEND 16 HZ ON THE left %%%%%%%%%%%%%%%   
                 [plotData_attend16hz_left,~,goodPos]=getDataSRC(c,2,1,2,3,folderName,folderLFP,occipital_list_right,useSingleElectrodeBadTrialList);
                 %%% baseline correction :
                  plotData_attend16hz_left_baselinecorrected=baselineCorrection(plotData_attend16hz_left,blPeriod);
                 %%spectrum estimation
                 [plotData_stimOn,plotData_baseline]=getDataForSpectrum(plotData_attend16hz_left_baselinecorrected,tERP,blPeriod);
                 [S1_attend16hz_left,mlogS1,dS1,t2,f2,average_S1_attend16hz_left] = getSpectrum(plotData_stimOn',timeVals,3,mtmParams,movingWin,0,blmin,blmax);
                 [S1_attend16hz_left_baseline,mlogS1,dS1,t2,f2,average_S1_attend16hz_left_baseline] = getSpectrum(plotData_baseline',timeVals,3,mtmParams,movingWin,0,blmin,blmax);
                 [S1,S1_power,no_of_trials]=getDifferenceSpectrum(S1_attend16hz_left,S1_attend16hz_left_baseline,average_S1_attend16hz_left,average_S1_attend16hz_left_baseline,useBaselineSpectrumFlag,findex,f2);
                 clear plotData_stimOn; clear plotData_baseline;
                 S1_attend16hz_left_spectrum{c}=S1;
                 S1_attend16hz_left_power(c)=S1_power;
                 S1_attend16hz_left_trials(c)=no_of_trials;
                 clear S1; clear S1_power; clear no_of_trials;
                 %%%%%%%%%%%%%%%%%%%%%%%%   UNATTEND 16 HZ ON THE LEFT %%%%%%%%%%%%%%% 
                  [plotData_unattend16hz_left,~,goodPos]=getDataSRC(c,3,1,1,3,folderName,folderLFP,occipital_list_right,useSingleElectrodeBadTrialList);
                 %%% baseline correction :
                  plotData_unattend16hz_left_baselinecorrected=baselineCorrection(plotData_unattend16hz_left,blPeriod);
                 %%spectrum estimation
                 [plotData_stimOn,plotData_baseline]=getDataForSpectrum(plotData_unattend16hz_left_baselinecorrected,tERP,blPeriod);
                 [S1_unattend16hz_left,mlogS1,dS1,t2,f2,average_S1_unattend16hz_left] = getSpectrum(plotData_stimOn',timeVals,3,mtmParams,movingWin,0,blmin,blmax);
                 [S1_unattend16hz_left_baseline,mlogS1,dS1,t2,f2,average_S1_unattend16hz_left_baseline] = getSpectrum(plotData_baseline',timeVals,3,mtmParams,movingWin,0,blmin,blmax);
                 [S1,S1_power,no_of_trials]=getDifferenceSpectrum(S1_unattend16hz_left,S1_unattend16hz_left_baseline,average_S1_unattend16hz_left,average_S1_unattend16hz_left_baseline,useBaselineSpectrumFlag,findex,f2);
                 clear plotData_stimOn; clear plotData_baseline;
                 S1_unattend16hz_left_spectrum{c}=S1;
                 S1_unattend16hz_left_power(c)=S1_power;
                 S1_unattend16hz_left_trials(c)=no_of_trials;
                 clear S1; clear S1_power; clear no_of_trials;
                 
                end 
                
            end
        end
    else %% second temporal frequency
        for j=1:2
            if j==1 % attend right 
                 for c=1:no_of_contrasts
                   %%%%%%%%%%%%%%%%%%%%%%%%   ATTEND 20 HZ ON THE RIGHT %%%%%%%%%%%%%%%   
                 [plotData_attend20hz_right,~,goodPos]=getDataSRC(c,3,1,1,3,folderName,folderLFP,occipital_list_left,useSingleElectrodeBadTrialList);
                 %%% baseline correction :
                  plotData_attend20hz_right_baselinecorrected=baselineCorrection(plotData_attend20hz_right,blPeriod);
                 %%spectrum estimation
                 [plotData_stimOn,plotData_baseline]=getDataForSpectrum(plotData_attend20hz_right_baselinecorrected,tERP,blPeriod);
                 [S1_attend20hz_right,mlogS1,dS1,t2,f2,average_S1_attend20hz_right] = getSpectrum(plotData_stimOn',timeVals,3,mtmParams,movingWin,0,blmin,blmax);
                 [S1_attend20hz_right_baseline,mlogS1,dS1,t2,f2,average_S1_attend20hz_right_baseline] = getSpectrum(plotData_baseline',timeVals,3,mtmParams,movingWin,0,blmin,blmax);
                 [S1,S1_power,no_of_trials]=getDifferenceSpectrum(S1_attend20hz_right,S1_attend20hz_right_baseline,average_S1_attend20hz_right,average_S1_attend20hz_right_baseline,useBaselineSpectrumFlag,findex,f2);
                 clear plotData_stimOn; clear plotData_baseline;
                 S1_attend20hz_right_spectrum{c}=S1;
                 S1_attend20hz_right_power(c)=S1_power;
                 S1_attend20hz_right_trials(c)=no_of_trials;
                 clear S1; clear S1_power; clear no_of_trials;
                 %%%%%%%%%%%%%%%%%%%%%%%%   UNATTEND 20 HZ ON THE RIGHT %%%%%%%%%%%%%%% 
                  [plotData_unattend20hz_right,~,goodPos]=getDataSRC(c,2,1,2,3,folderName,folderLFP,occipital_list_left,useSingleElectrodeBadTrialList);
                 %%% baseline correction :
                  plotData_unattend20hz_right_baselinecorrected=baselineCorrection(plotData_unattend20hz_right,blPeriod);
                 %%spectrum estimation
                 [plotData_stimOn,plotData_baseline]=getDataForSpectrum(plotData_unattend20hz_right_baselinecorrected,tERP,blPeriod);
                 [S1_unattend20hz_right,mlogS1,dS1,t2,f2,average_S1_unattend20hz_right] = getSpectrum(plotData_stimOn',timeVals,3,mtmParams,movingWin,0,blmin,blmax);
                 [S1_unattend20hz_right_baseline,mlogS1,dS1,t2,f2,average_S1_unattend20hz_right_baseline] = getSpectrum(plotData_baseline',timeVals,3,mtmParams,movingWin,0,blmin,blmax);
                 [S1,S1_power,no_of_trials]=getDifferenceSpectrum(S1_unattend20hz_right,S1_unattend20hz_right_baseline,average_S1_unattend20hz_right,average_S1_unattend20hz_right_baseline,useBaselineSpectrumFlag,findex,f2);
                 clear plotData_stimOn; clear plotData_baseline;
                 S1_unattend20hz_right_spectrum{c}=S1;
                 S1_unattend20hz_right_power(c)=S1_power;
                 S1_unattend20hz_right_trials(c)=no_of_trials;
                 clear S1; clear S1_power; clear no_of_trials;  
                 end    
            else   
                 for c=1:no_of_contrasts    
                    %%%%%%%%%%%%%%%%%%%%%%%%   ATTEND 20 HZ ON THE LEFT %%%%%%%%%%%%%%%   
                 [plotData_attend20hz_left,~,goodPos]=getDataSRC(c,3,1,2,3,folderName,folderLFP,occipital_list_right,useSingleElectrodeBadTrialList);
                 %%% baseline correction :
                  plotData_attend20hz_left_baselinecorrected=baselineCorrection(plotData_attend20hz_left,blPeriod);
                 %%spectrum estimation
                 [plotData_stimOn,plotData_baseline]=getDataForSpectrum(plotData_attend20hz_left_baselinecorrected,tERP,blPeriod);
                 [S1_attend20hz_left,mlogS1,dS1,t2,f2,average_S1_attend20hz_left] = getSpectrum(plotData_stimOn',timeVals,3,mtmParams,movingWin,0,blmin,blmax);
                 [S1_attend20hz_left_baseline,mlogS1,dS1,t2,f2,average_S1_attend20hz_left_baseline] = getSpectrum(plotData_baseline',timeVals,3,mtmParams,movingWin,0,blmin,blmax);
                 [S1,S1_power,no_of_trials]=getDifferenceSpectrum(S1_attend20hz_left,S1_attend20hz_left_baseline,average_S1_attend20hz_left,average_S1_attend20hz_left_baseline,useBaselineSpectrumFlag,findex,f2);
                 clear plotData_stimOn; clear plotData_baseline;
                 S1_attend20hz_left_spectrum{c}=S1;
                 S1_attend20hz_left_power(c)=S1_power;
                 S1_attend20hz_left_trials(c)=no_of_trials;
                 clear S1; clear S1_power; clear no_of_trials;
                 %%%%%%%%%%%%%%%%%%%%%%%%   UNATTEND 20 HZ ON THE LEFT %%%%%%%%%%%%%%% 
                  [plotData_unattend20hz_left,~,goodPos]=getDataSRC(c,2,1,1,3,folderName,folderLFP,occipital_list_left,useSingleElectrodeBadTrialList);
                 %%% baseline correction :
                  plotData_unattend20hz_left_baselinecorrected=baselineCorrection(plotData_unattend20hz_left,blPeriod);
                 %%spectrum estimation
                 [plotData_stimOn,plotData_baseline]=getDataForSpectrum(plotData_unattend20hz_left_baselinecorrected,tERP,blPeriod);
                 [S1_unattend20hz_left,mlogS1,dS1,t2,f2,average_S1_unattend20hz_left] = getSpectrum(plotData_stimOn',timeVals,3,mtmParams,movingWin,0,blmin,blmax);
                 [S1_unattend20hz_left_baseline,mlogS1,dS1,t2,f2,average_S1_unattend20hz_left_baseline] = getSpectrum(plotData_baseline',timeVals,3,mtmParams,movingWin,0,blmin,blmax);
                 [S1,S1_power,no_of_trials]=getDifferenceSpectrum(S1_unattend20hz_left,S1_unattend20hz_left_baseline,average_S1_unattend20hz_left,average_S1_unattend20hz_left_baseline,useBaselineSpectrumFlag,findex,f2);
                 clear plotData_stimOn; clear plotData_baseline;
                 S1_unattend20hz_left_spectrum{c}=S1;
                 S1_unattend20hz_left_power(c)=S1_power;
                 S1_unattend20hz_left_trials(c)=no_of_trials;
                 clear S1; clear S1_power; clear no_of_trials; 
                    
                 end 
                
            end 
        end 
    end
end
 
%% FIGURES :
%%%% attend 16 vs unattend 16 graphs :

plotHandles_barPlot_F1=getPlotSSVEPHandles(2,no_of_contrasts);
cMidindex=round(max(no_of_contrasts)/2);
%%Contrast Labels:
cIndexValsUnique=cindex;
for c=1:length(cIndexValsUnique)
    clabel{c}=num2str([cIndexValsUnique(c)]);
end
   clabel{length(cIndexValsUnique)+1}='all';

for c=1:no_of_contrasts
subplot(plotHandles_barPlot_F1(1,c));
y=[S1_attend16hz_right_power(c) S1_unattend16hz_right_power(c)];
bar(1,y(1),'FaceColor','b');hold on;bar(2,y(2),'FaceColor','r');ylim([-3 3]);
text(0.1,0.95,['n = ' num2str(S1_attend16hz_right_trials(c))],'unit','normalized','fontsize',9);
text(0.55,0.95,['n = ' num2str(S1_unattend16hz_right_trials(c))],'unit','normalized','fontsize',9);
if c==cMidindex
    title( 'Attend 16 hz right vs unattend 16 hz right');
end
if c==no_of_contrasts
      legend('Attend','UnAttend','Location','best');
end
if c==1
    ylabel('SSVEP Amplitude at 16 Hz');
end
subplot(plotHandles_barPlot_F1(2,c));
y=[S1_attend16hz_left_power(c) S1_unattend16hz_left_power(c)];
bar(1,y(1),'FaceColor','b');hold on;bar(2,y(2),'FaceColor','r');ylim([-3 3]);xlabel([clabel{c}]);
text(0.1,0.95,['n = ' num2str(S1_attend16hz_left_trials(c))],'unit','normalized','fontsize',9);
text(0.55,0.95,['n = ' num2str(S1_unattend16hz_left_trials(c))],'unit','normalized','fontsize',9);
 if c==cMidindex
     title( 'Attend 16 hz left vs unattend 16 hz left');
 end 
end
% %% fft for the above
figure;
plotHandles_FFTPlot_F1=getPlotSSVEPHandles(2,no_of_contrasts);
for c=1:no_of_contrasts
subplot(plotHandles_FFTPlot_F1(1,c));
plot(f2,S1_attend16hz_right_spectrum{c},'b');
hold on;
plot(f2,S1_unattend16hz_right_spectrum{c},'r');
if c==cMidindex
    title( 'Attend 16 hz right vs unattend 16 hz right');
end
if c==no_of_contrasts
      legend('Attend','UnAttend','Location','best');
end
if c==1
    ylabel('SSVEP Amplitude at 16 Hz');
end
xlim([0 200]);
ylim([0 60]);
subplot(plotHandles_FFTPlot_F1(2,c));plot(f2,S1_attend16hz_left_spectrum{c},'b');
hold on;
plot(f2,S1_unattend16hz_left_spectrum{c},'r');
if c==cMidindex
     title( 'Attend 16 hz left vs unattend 16 hz left');
 end 
xlim([0 200]);
ylim([0 60]);
end

figure;
plotHandles_barPlot_F2=getPlotSSVEPHandles(2,no_of_contrasts);
for c=1:no_of_contrasts
subplot(plotHandles_barPlot_F2(1,c));
y=[S1_attend20hz_right_power(c) S1_unattend20hz_right_power(c)];
bar(1,y(1),'FaceColor','b');hold on;bar(2,y(2),'FaceColor','r');ylim([-3 3]);
text(0.1,0.95,['n = ' num2str(S1_attend20hz_right_trials(c))],'unit','normalized','fontsize',9);
text(0.55,0.95,['n = ' num2str(S1_unattend20hz_right_trials(c))],'unit','normalized','fontsize',9);
if c==cMidindex
    title( 'Attend 20 hz right vs unattend 20 hz right');
end
if c==no_of_contrasts
      legend('Attend','UnAttend','Location','best');
end
if c==1
    ylabel('SSVEP Amplitude at 20 Hz');
end
subplot(plotHandles_barPlot_F2(2,c));
y=[S1_attend20hz_left_power(c) S1_unattend20hz_left_power(c)];
bar(1,y(1),'FaceColor','b');hold on;bar(2,y(2),'FaceColor','r');ylim([-3 3]);xlabel([clabel{c}]);
text(0.1,0.95,['n = ' num2str(S1_attend20hz_left_trials(c))],'unit','normalized','fontsize',9);
text(0.55,0.95,['n = ' num2str(S1_unattend20hz_left_trials(c))],'unit','normalized','fontsize',9);
 if c==cMidindex
     title( 'Attend 20 hz left vs unattend 20 hz left');
 end 
end

% %% fft for the above 
figure;
plotHandles_FFTPlot_F2=getPlotSSVEPHandles(2,no_of_contrasts); 
for c=1:no_of_contrasts
subplot(plotHandles_FFTPlot_F2(1,c));
plot(f2,S1_attend20hz_right_spectrum{c},'b');
hold on;
plot(f2,S1_unattend20hz_right_spectrum{c},'r');
if c==cMidindex
    title( 'Attend 20 hz right vs unattend 16 hz right');
end
if c==no_of_contrasts
      legend('Attend','UnAttend','Location','best');
end
if c==1
    ylabel('SSVEP Amplitude at 20 Hz');
end
xlim([0 200]);
ylim([0 60]);
subplot(plotHandles_FFTPlot_F2(2,c));plot(f2,S1_attend20hz_left_spectrum{c},'b');
hold on;
plot(f2,S1_unattend20hz_left_spectrum{c},'r');
if c==cMidindex
     title( 'Attend 20 hz left vs unattend 20 hz left');
 end 
xlim([0 200]);
ylim([0 60]);
end           

end






function [Data,trialNums,goodTrials]=getDataSRC(c,t,e,a,s,folderName,folderLFP,occipital_list,useSingleElectrodeBadTrialList)

folderExtract = fullfile(folderName,'extractedData');
    folderSegment = fullfile(folderName,'segmentedData');
    
    [parameterCombinations] = loadParameterCombinations(folderExtract);
    try
        load(fullfile(folderSegment,'badTrials.mat'));
    catch
        disp('No bad trials')
        badTrials = [];
    end
    trialNums = cell2mat(parameterCombinations(c,t,e,a,s));
    
% Extraction
     
%         for ele=1:length(occipital_list)
            if useSingleElectrodeBadTrialList
    
            goodTrials = setdiff(trialNums,allBadTrials{occipital_list});
     else 
         
         goodTrials = setdiff(trialNums,badTrials);
            end 
        analogData=loadAnalogData(fullfile(folderLFP,['elec' num2str(occipital_list) '.mat']));
%          Data(ele,:,:)=analogData(goodTrials,:);
       
%         end
%         
%         Data=mean(Data);
%         Data=squeeze(Data);
       analogData=loadAnalogData(fullfile(folderLFP,['elec' num2str(occipital_list) '.mat']));
       Data=analogData(goodTrials,:);
    
     
end 

function stimResults = loadStimResults(folderExtract)
load (fullfile(folderExtract,'stimResults'));
end


function plotData_baselinecorrected=baselineCorrection(plotData_uncorrected,blPeriod)
if ~exist('blPeriod','var')
    blPeriod=0;
end
for ei=1:size(plotData_uncorrected,1)
         baseline_factor=mean(plotData_uncorrected(ei,blPeriod));
         plotData_baselinecorrected(ei,:)=plotData_uncorrected(ei,:)-baseline_factor;
end 
end

function [plotData_stimOn,plotData_baseline]=getDataForSpectrum(plotData,tERP,blPeriod)

if ~exist('tERP','var') tERP=0; end
if ~exist('blPeriod','var') blPeriod=0; end

plotData_stimOn=plotData(:,tERP);
plotData_baseline=plotData(:,blPeriod);
end



function [S1,S1_power,no_of_trials]=getDifferenceSpectrum(S1_stimOn,S1_baseline,averageS1,averageS1_baseline,useBaselineSpectrumFlag,findex,f2)

if ~exist('useTrialWiseSubtraction','var') 
    useTrialWiseSubtraction=1; 
end
if ~exist('useBaselineSpectrumFlag','var') 
    useBaselineSpectrumFlag=0;
end
 
if useTrialWiseSubtraction && useBaselineSpectrumFlag
    
    %% baseline spectrum is calculated for each trial and subtracted from spectrum of corresponding trial and
    %% averaged after that. 
    
 for i=1:size(S1_stimOn,2)
      S1(:,i)=S1_stimOn(:,i)-S1_baseline(:,i);
      S1_power(i)=S1((findex(2)==f2),i);
      
 end
    S1=mean(S1,2);
    S1_power=mean(S1_power);
    
else if (useTrialWiseSubtraction==0) && useBaselineSpectrumFlag
        %% here baseline spectrum is averaged across trials and subtracted from each trial ( all the trials are subjected to
        %% the same baseline spectrum 
        
        S1=averageS1-averageS1_baseline;
        S1_power=averageS1(findex==f2)-averageS1_baseline(findex==f2);
       
            
    else 
            S1=S1_stimOn;
            S1_power=S1_stimOn(findex==f2);
    end    
end   
no_of_trials=size(S1_stimOn,2);
    

end 

