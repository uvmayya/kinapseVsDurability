function [samplingEff,meanOverlap] = getSamplingAndOverlapStats_v2(datacellInit)
% calculates avg. sampling efficiency and overlap time-course for the cells that were attached on spots and didn't have another attached cell on the same spot. 
%   Detailed explanation goes here

% copied this part from getTimeBlockMotilityData.m in 041414/anals/072415_motilityAnals
%frameInterval=28.94; % in seconds
expName='070617_2o2co_088_1_sampOverlapPerTrack.mat';
flurThresh=50; %used to filter out attached portions of tracks

% to generate the data-specific parameters needed for getting DIC outlines and sampling efficiency
dirstring = '2o2co_088_1/';
%dirstring = '../nve_allCh/';
DICstartImg = 2; 
spotImg = imread('spot_2o2co_088_1.tif'); %[dirstring,'mem_4minOn_t001_c002.tif']
%spotImg = imread([dirstring,'nve_13minOn_t001_c002.tif']); 
cyclesize = 2;
numimgs = 180; %480; %100
videocell = imgfolder2videocell(dirstring, DICstartImg, cyclesize, numimgs);
halfCropSize = 23; % Initially I was using 22 for nve and 23 for mem, but then decided to keep the same (i.e. 23) for both
detectParams = [1.0, 0.1, 5, 13, 10, 15, 5, 1]; % need to brighten DIC for mem
%detectParams = [1.1, 0.10, 4, 14, 10, 15, 6, 0];
% detection parameters explained: 
umPerPixel = 0.414; 

minLength=30; 
proxThresh=13; % in um
datacell=selectArrestedPortionsOnSpots(datacellInit,flurThresh,minLength,proxThresh); 

samplingEff=zeros(length(datacell),1); 
meanOverlap=zeros(length(datacell),3);
for cp=1:length(datacell)
    trkLength=length(datacell{cp}(:,1));
    nChunks=floor(trkLength/20); % considering tracks in blocks of 20 or 40 steps
    sampEff=zeros(nChunks,1); overlap=zeros(nChunks,3);
    for chunk=1:nChunks
        startPos=((chunk-1)*20)+1;
        endPos=chunk*20;
        trackChunk=datacell{cp}(startPos:endPos,:);
        stFrame = datacell{cp}(1, 1) + startPos -1;
        [sampEff(chunk),overlap(chunk,1:3)] = getSampEffOverlapForTrack(trackChunk, videocell, spotImg, stFrame, halfCropSize, detectParams, umPerPixel);
    end    
    samplingEff(cp)=mean(sampEff);
    meanOverlap(cp,1)=mean(overlap(:,1)); % fractional overlap
    meanOverlap(cp,2)=mean(overlap(:,2)); % overlap relative to spot area
    meanOverlap(cp,3)=mean(overlap(:,3)); % cell area relative to spot area
    fprintf('Sampling efficiency and overlap calculated for cell-track: %d/%d \n',cp,size(datacell,2));
end


save(expName,'samplingEff','meanOverlap');