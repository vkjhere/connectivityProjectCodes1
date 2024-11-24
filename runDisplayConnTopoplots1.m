% plotFigures
comparisonStr = 'paired';
protocolName = 'G1';
analysisChoice = 'st';

badEyeCondition = 'ep';
badTrialVersion = 'v8';

freqRangeList{1} = [8 12];
%  freqRangeList{2} = [24 34];
freqRangeList{2} = [20 25];
% freqRangeList{3} = [30 80];
freqRangeList{3} = [60  98];

electrodeList = [16 17 18 48]; %(O1-Oz-O2-POz);  [14 44 47];%(P3-P1-PO3)     % [19 49 52];%(P4-PO4-P2) ;           % ; % 
connMethod = 'ppc';
displayDataFlag = 1;

axisRangeList{1} = [0 1]; axisRangeName{1} = 'YLims';
axisRangeList{2} = [0 1]; axisRangeName{2} = 'cLims (topo)';

cutoffList = [2 30];
useMedianFlag = 0;

if strcmp(comparisonStr,'paired')
    pairedSubjectNameList = getPairedSubjectsBK1;
    subjectNameLists{1} = pairedSubjectNameList(:,1);
    subjectNameLists{2} = pairedSubjectNameList(:,2);
%  subjectNameLists{1} = {'013AR'};
%  subjectNameLists{2} = {'064PK'};
    pairedDataFlag      = 1;
else
    [~, meditatorList, controlList] = getGoodSubjectsBK1;
    subjectNameLists{1} = meditatorList;
    subjectNameLists{2} = controlList;
    pairedDataFlag      = 0;
end
hAllPlots = [];


[connDataToReturn,goodSubjectNameListsToReturn,topoplotDataToReturn,freqVals] = displayConnTopoplotsAllSubjects1(subjectNameLists,protocolName,analysisChoice,electrodeList,connMethod,badEyeCondition,badTrialVersion,freqRangeList,axisRangeList,cutoffList,useMedianFlag,hAllPlots,pairedDataFlag,1);
