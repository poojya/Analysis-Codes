% run finBadTrialsEEG
% Vinay Shirhatti, 13 Oct 2014
%==========================================================================



clear all; clc;
% Choose the protocols (the indices correspond to those in listProtocols.m)
extractTheseIndices = 79;
subjectName = 'Human'; gridType = 'EEG';

% Get the details for protocols
[subjectNames,expDates,protocolNames,stimTypes] = eval(['allProtocols' upper(subjectName(1)) subjectName(2:end) gridType]);

% Check the OS and set paths accordingly
if ispc
    folderSourceString = 'W:\';
else
    folderSourceString = '/media/store/';
end

% define grid
 gridType = 'EEG';
% gridType = 'Microelectrode';




% occipitalElec = [9,10,45,46,59,60,63,64];
% parietalElec = [7,8,15,16,19,37,38,51,52];

% occipitalElec = [20,21,22,25,26,27];
% parietalElec = [7,8,15,16,19,37,38,51,52];

% occipitalElec = [1:3];

% Main loop to determine the bad trials
for i = 1:length(extractTheseIndices)
    index = extractTheseIndices(i);
    disp(['Find bad trials for index: ' num2str(index)]);
    subjectName = subjectNames{index};
    expDate = expDates{index};
    protocolName = protocolNames{index};
    
    % check these electrodes to decide the overall bad trials
%     checkTheseElectrodes = [occipitalElec,parietalElec];
%     checkTheseElectrodes = occipitalElec;
    checkTheseElectrodes = [29 30 31];
%     checkTheseElectrodes = 1:96;
    
    % Set limits/criteria
    threshold = 6;
    maxLimit = 100; minLimit = -100;
    saveData = 1;
    showTrials = 0;
    showElectrodes = 43;
    checkPeriod = [-0.4 0.6];
    
    [~,~,~] = findBadTrialsEEG(subjectName,expDate,...
        protocolName,folderSourceString,gridType,...
        checkTheseElectrodes,threshold,maxLimit,minLimit,saveData,showTrials,checkPeriod);
%     
%     [~,~] = findBadTrialsWithLFPv2(subjectName,expDate,...
%         protocolName,folderSourceString,gridType,...
%         checkTheseElectrodes,threshold,maxLimit,showElectrodes,minLimit,saveData,checkPeriod);
end

