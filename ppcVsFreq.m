% To plot PPC vs freq and cos(theta) of different elecs wrt seeds

comparisonStr = 'paired';
protocolName = 'G1';
analysisChoice = 'st';

badEyeCondition = 'ep';
badTrialVersion = 'v8';

electrodeList =  [16 17 18 48];%(O1-Oz-O2-POz); [14 44 47];%(P3-P1-PO3) %    ; [19 49 52];%(P4-PO4-P2);           % ; % 
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
%     subjectNameLists{1} = {'013AR'}; %star meditator
%     subjectNameLists{2} = {'064PK'}; %paired control

    pairedDataFlag      = 1;
else
    [~, meditatorList, controlList] = getGoodSubjectsBK1;
    subjectNameLists{1} = meditatorList;
    subjectNameLists{2} = controlList;
    pairedDataFlag      = 0;
end
hAllPlots = [];

if ~exist('protocolName','var');          protocolName='G1';            end
if ~exist('analysisChoice','var');        analysisChoice='st';          end
if ~exist('electrodeList','var');         electrodeList = [];           end
if ~exist('connMethod','var');            connMethod = 'ppc';           end
if ~exist('badEyeCondition','var');       badEyeCondition='ep';         end
if ~exist('badTrialVersion','var');       badTrialVersion='v8';         end
if ~exist('freqRangeList','var')
    freqRangeList{1} = [8 13]; % alpha
    freqRangeList{2} = [20 30]; % SG 
    freqRangeList{3} = [35 65]; % FG
end
if ~exist('axisRangeList','var')
    axisRangeList{1} = [0 100];
    axisRangeList{2} = [-2.5 2.5];
    axisRangeList{3} = [-1.5 1.5];
end
if ~exist('cutoffList','var')
    cutoffList = [3 30];
end

if ~exist('useMedianFlag','var');         useMedianFlag = 0;            end
if ~exist('hAllPlots','var');             hAllPlots = [];               end
if ~exist('pairedDataFlag','var');        pairedDataFlag = 0;           end
if ~exist('displayDataFlag','var');       displayDataFlag = 1;          end

numFreqRanges = length(freqRangeList);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Display options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
displaySettings.fontSizeLarge = 10;
displaySettings.tickLengthMedium = [0.025 0];
displaySettings.colorNames(1,:) = [0.8 0 0.8];      % Purple 
displaySettings.colorNames(2,:) = [0.25 0.41 0.88]; % Cyan
titleStr{1} = 'Meditators';
titleStr{2} = 'Controls';

cLimsTopo = axisRangeList{2};
binRange = [-0.5 0.5]; % why?
%%%%%%%%%%%%%%%%%%%%%%%%%% Get topoplot info %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
capType = 'actiCap64_UOL';
x = load([capType '.mat']); 
montageChanlocs = x.chanlocs; % channel locations x,y,z etc.

saveFolderName = 'savedData1';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[connData,freqVals,connDataBinwise,binnedCenters] = getConnDataAllSubjectsV(subjectNameLists,electrodeList,connMethod,badEyeCondition,badTrialVersion,protocolName,analysisChoice,cutoffList,pairedDataFlag,saveFolderName,montageChanlocs);


%%%%%%%%%%%%%%%%%%%%%%% plots: PPC vs freq and cos(theta)%%%%%%%%%%%%%%%%%%%%%
% to define common scale of coorbar for the first 2 plots. Difference
% plot's color scale is kept free.
meanData1 = squeeze(mean((connDataBinwise{1}),1, 'omitnan')); % meditators 
meanData2 = squeeze(mean((connDataBinwise{2}),1, 'omitnan')); % controls
diffData = squeeze(mean((connDataBinwise{1} - connDataBinwise{2}),1, 'omitnan'));
% Determine the common limits for colorbar scale
cmin = min([min(meanData1(:)), min(meanData2(:))]);
cmax = max([max(meanData1(:)), max(meanData2(:))]);



figure(1)
sgtitle(['PPC wrt elecs ', num2str(electrodeList) , ' during ',protocolName]);

n = length(freqVals)/20; % Change 'n' based on how many labels you want
uniform_xticks = freqVals(1:n:end);

% plotting ppc just for nearby electrodes (back hemisphere) : cos(dist) :
% [0 1]
binnedCenters = binnedCenters(5:8); 
connDataBinwise{1} = connDataBinwise{1}(:,5:8,:);
connDataBinwise{2} = connDataBinwise{2}(:,5:8,:);

subplot(311) % meditators
pcolor(freqVals, binnedCenters, squeeze(mean((connDataBinwise{1}),1, 'omitnan')));
colormap("jet"); 
set(gca, 'XScale', 'log');
set(gca, 'XTick', uniform_xticks);
set(gca, 'XTickLabel', arrayfun(@num2str, uniform_xticks, 'UniformOutput', false));
shading interp;  axis tight;   clim([cmin cmax]);
%xlim([20 40]);
xlabel('frequency'); ylabel('cos (theta)'); title('Meditators');

subplot(312) %control
pcolor(freqVals, binnedCenters, squeeze(mean((connDataBinwise{2}),1, 'omitnan')));
shading interp;  axis tight;   clim([cmin cmax]);
set(gca, 'XScale', 'log');
set(gca, 'XTick', uniform_xticks);
set(gca, 'XTickLabel', arrayfun(@num2str, uniform_xticks, 'UniformOutput', false));
%xlim([ 20 40]);
xlabel('frequency'); ylabel('cos (theta)'); title('Controls');

subplot(313) % diff
pcolor(freqVals, binnedCenters, squeeze(mean((connDataBinwise{1}-connDataBinwise{2}),1, 'omitnan'))); 
shading interp; colorbar; axis tight; clim([-0.04 0.04]);
set(gca, 'XScale', 'log');
set(gca, 'XTick', uniform_xticks);
set(gca, 'XTickLabel', arrayfun(@num2str, uniform_xticks, 'UniformOutput', false));
%xlim([ 20 45]);
xlabel('frequency'); ylabel('cos (theta)'); title('diff.');

%%%%%%%%%%%%%%%%%%%%% functions %%%%%%%%%%%%%%%%%%%%%%%%%%
function [connData,freqVals,connDataBinwise,binnedCenters] = getConnDataAllSubjectsV(subjectNameLists,electrodeList,connMethod,badEyeCondition,badTrialVersion,protocolName,analysisChoice,cutoffList,pairedDataFlag,saveFolderName,montageChanlocs)

%%%%% Discretize connectivty into bins depending on distance from seed %%%%
binWidth = 0.25;
binEdges = -1:binWidth:1;
nbins = length(binEdges)-1;
binnedCenters = binEdges(1:end-1)+(binWidth/2);
loc = getElecLocAngles(montageChanlocs);

numElectrodes = length(electrodeList); % = 3 or 4 usually
binnedIndicesAllElectrodes = cell(numElectrodes,nbins);

for e=1:numElectrodes % 3 or 4 usually
    dist = sqrt((angl_dist(loc.azi(electrodeList(e)),loc.azi,'a')).^2+(angl_dist(loc.ele(electrodeList(e)),loc.ele,'e')).^2);
    fitx = cos((dist/180)*pi);
    binned_fitx = discretize(fitx,binEdges);
    
    for b = 1:nbins % number of bins
        binnedIndicesAllElectrodes{e,b} = find(binned_fitx == b); % w r t each of the given 3 or 4 elecs, 
        % it calculates the electrodes which have cos(dist) in different bins like (0.5 0.75)
        % bin_index = 8 means, cos (dist) = 1; so the 2 electrodes are close to
        % each other.
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize
badSubjectList = cell(1,2);
connDataTMP = cell(1,2);
connDataBinwiseTMP = cell(1,2);

for i=1:2
    for j=1:length(subjectNameLists{i}) % 2 lists: med and con
        subjectName = subjectNameLists{i}{j};

        tmpData = load(fullfile(saveFolderName,subjectName,[protocolName '_' badEyeCondition '_' badTrialVersion '_' connMethod]));
        numGoodTrials = tmpData.numGoodTrials;

        if numGoodTrials<cutoffList(2)
            badSubjectList{i}(j) = 1;
        else
            if strcmp(analysisChoice,'bl')
                connDataTMP2 = tmpData.connPre(electrodeList,:,:); % conn 64*64 elecs for the given subject
                freqVals = tmpData.freqPre;
            elseif strcmp(analysisChoice,'st')
                connDataTMP2 = tmpData.connPost(electrodeList,:,:);
                freqVals = tmpData.freqPost;
            else
                connDataTMP2 = (tmpData.connPre(electrodeList,:,:) + tmpData.connPost(electrodeList,:,:))/2;
                freqVals = tmpData.freqPost;
            end

            numGoodElectodes = trace(~isnan(squeeze(connDataTMP2(:,electrodeList,1))));

            if numGoodElectodes >= cutoffList(1)
                badSubjectList{i}(j) = 0;
                connDataTMP{i}{j} = squeeze(mean(connDataTMP2,1,'omitnan'));
                
                connDataBinwiseTMP2 = zeros(numElectrodes,nbins,length(freqVals));
                for e=1:numElectrodes
                    for b=1:nbins
                         connDataBinwiseTMP2(e,b,:) = squeeze(mean(connDataTMP2(e,binnedIndicesAllElectrodes{e,b},:),2,'omitnan'));
                    end
                end
                connDataBinwiseTMP{i}{j} = squeeze(mean(connDataBinwiseTMP2,1,'omitnan'));
            else
                badSubjectList{i}(j) = 1;
            end
        end
    end
end

% Remove bad subjects
connData = cell(1,2);
connDataBinwise = cell(1,2);

for i=1:2
    if pairedDataFlag
        badSubjectPos = find(sum(cell2mat(badSubjectList')));
    else
        badSubjectPos = find(badSubjectList{i});
    end
    x1 = connDataTMP{i};
    x1(badSubjectPos)=[];
    x2 = connDataBinwiseTMP{i};
    x2(badSubjectPos)=[];

    numSubjects = length(x1);
    for j=1:numSubjects
        connData{i}(j,:,:) = x1{j};
        connDataBinwise{i}(j,:,:) = x2{j};
    end
end
end

%%%%%%%%%%%% functions %%%%%%%%%%%%%%%%%%%

function out_theta = angl_dist(in_theta_ref,in_theta,val)
if(strcmp(val,'a')) % azimuth (addressing Cz issue)
    if(in_theta_ref > 90)
    in_theta(14) = 90;
    elseif(in_theta_ref < -90)
        in_theta(14) = -90;
    end
end
in_theta = abs(in_theta-in_theta_ref);
out_theta = zeros(1,length(in_theta));
for i=1:length(in_theta)
    if(in_theta(i) > 180)
        out_theta(i) = 360 - in_theta(i); % to get shortest angular distance
    else
        out_theta(i) = in_theta(i);
    end
end
if(strcmp(val,'a')) % azimuth (addressing Cz issue)
    out_theta(14) = 0;
end
end