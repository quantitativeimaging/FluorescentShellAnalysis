% A script to track ellipsoidal spores by post-processing ELM output data
% 
% 24 Nov 2015
%
% Version: 1. This is a first attempt and has not yet been validated.
% 
% Method: Find nearest neighbour in previous frame, and link to this frame
%         Only accept if a single candidate exists (not zero or many) 
%           within an acceptable search radius. 
%         Only keep candidates where a single option exists in all frames
% 

% 1. SETUP.
%    Which folder of .mat files to read?
%    Which folder to write frames of video output?
%    Max jump distance to acceptable nearest neighbours

inDir   = 'C:\Users\user\Documents\Projects\2015_IIB_Spores\Data\19_11_15_results\mats\';
file1   = '\SleL_GFP_2a.mat';

outDir = 'C:\Users\user\Documents\Projects\2015_IIB_Spores\Data\19_11_15_results\vid\';

jumpRadSq = 4^2; % Max unambiguous jump distance. 

listOfFiles   = dir([inDir,'*.mat']); % Get all .mat files
numberOfFiles = length(listOfFiles);


% 2. ANALYSIS 
%
% Load the first frame of data and find out how many spores there are
load([inDir,file1]);
numberSporesFrame1 = length(listFittedRow);
listLength = numberSporesFrame1 + 10;

listTrackedIndex = -ones(listLength,numberOfFiles);
listTrackedIndex(1:numberSporesFrame1, 1) = listFittedInd;

% Get the position and index numbers of spores in first frame
listRowOld = listFittedRow;
listColOld = listFittedCol;
listIndOld = listFittedInd;

% Find the tracked number associated with spores in later frames
% Set older entries to -1 if not found in subsequent frame
 for lpFile = 2:numberOfFiles
    lpFile % Show progress
     
    load([inDir,listOfFiles(lpFile).name]);
    listRowNew = listFittedRow;
    listColNew = listFittedCol;
    listIndNew = listFittedInd;
    
    for lpCand = 1:length(listRowNew)
       dispsSq = (listRowNew(lpCand) - listRowOld).^2 + ...
                 (listColNew(lpCand) - listColOld).^2;
       
       posCands= dispsSq < jumpRadSq;
       
       if(sum(posCands)==1);
           [m,Indx] = min(dispsSq); % Find closest previous spore
           candOldInd = listTrackedIndex(Indx, lpFile-1);
           listTrackedIndex(lpCand, lpFile) = candOldInd;
       end
    end
    
    % Store current positions and indexes for 
    listRowOld = listRowNew;
    listColOld = listColNew;
    listIndOld = listIndNew;
 end
 
 
%% Show tags on recon^s. Only number spores found in all frames
flagVidRecon =1;
figure(1)
for lpFile = 1:numberOfFiles
   load([inDir,listOfFiles(lpFile).name]);
   
   if(flagVidRecon)
    sim_all; % DANGER: Deletes some listFittedXXX for quality control 
    figure(1)
    imagesc(imSim'); % Note imagesc (re)-defines y-axis to be down-pointing
    colormap(gray)
    load([inDir,listOfFiles(lpFile).name]); % RELOAD overwrites by sim_all
   else 
    scatter(listFittedCol*scaleFactor, listFittedRow*scaleFactor, '+g')
   end
   hold on
     scatter(listFittedCol*scaleFactor, listFittedRow*scaleFactor, '+g')
 
     b = listTrackedIndex(1:length(listFittedRow),lpFile);
     c = num2str(b);
     c((b==-1), :) = ' '; % Don't show -1 errors.
     dx = 5; 
     dy = -1; % displacement so the text does not overlay the data point
     text((listFittedCol+dx)*scaleFactor,(listFittedRow+dy)*scaleFactor,...
          c, 'color', 'g', 'fontSize', 10);
   hold off
     
   pause(2)
   
   if(1)% select here whether to write some output
       set(1,'Position',[100,100,800,600]); % 800 px wide, 600 high
       xlim([1001 1800]); % Choose area of interest
       ylim([[1 600]]);
       F   = getframe(1);
       imF = F.cdata;
       imwrite(imF, [outDir, 'vid', int2str(lpFile), '.png'])
   end
end

%% Plot graph of axis length changes
%
% Inputs: listTrackedIndex and .Mats of listFittedRad, listFittedAR
% Method: collate list of spore properties with a double for loop

sporeInd = listTrackedIndex(1,:);
sporeInd = sporeInd(sporeInd ~= -1); % Remove surplus -1 padding

sporeAxA = -ones(length(sporeInd), numberOfFiles); % Semi-minor Axis Length
sporeAxB = -ones(length(sporeInd), numberOfFiles); % Semi-major Axis Length

for lpFile = 1:numberOfFiles
  load([inDir,listOfFiles(lpFile).name]);
  for lpSpore = 1:length(listFittedRad) 
      if( listTrackedIndex(lpSpore,lpFile)~=(-1) ) % Don't assign if can't
          %lpSpore
          %listTrackedIndex(lpSpore,lpFile)~=(-1)
      sporeAxA(listTrackedIndex(lpSpore,lpFile), lpFile) = listFittedRad(lpSpore);
      sporeAxB(listTrackedIndex(lpSpore,lpFile), lpFile) = listFittedRad(lpSpore).*listFittedAR(lpSpore);
      end
  end
end

% Note that zero seems to be the "unassigned" code in listFittedRad
tData = 1:numberOfFiles;
figure(2)
plot(tData(1,sporeAxA(1,:)~=0), sporeAxA(1,sporeAxA(1,:)~=0));
hold on
plot(tData(1,sporeAxB(1,:)~=0), sporeAxB(1,sporeAxA(1,:)~=0), 'r--');
for lpSpore = 1:size(sporeAxA,1) % for each row of data
     plot(-5+5*tData(sporeAxA(lpSpore,:)~=0), ...
          74*sporeAxA(lpSpore,sporeAxA(lpSpore,:)~=0), 'b')
     plot(-5+5*tData(sporeAxB(lpSpore,:)~=0), ...
          74*sporeAxB(lpSpore,sporeAxB(lpSpore,:)~=0), 'r--')
end
hold off

legend('Semi-minor axis','Semi major axis')
xlabel('time, min')
ylabel('length, nm')
title('Non-germinating spores')
set(gcf,'color','w')
