% Image analysis of spores
%   Find model parameters for spheres or ellipsoids to fit image data
%   
% Eric Rees March 2015. 
%   Exact solution for spherical shell: David Rees.
%   Other models and fitting methods: Eric Rees
%   Example spore image data: Graham Christie
%
% Licence: LGPL Version 3+, http://www.gnu.org/licenses/lgpl.html 
%
% Website: http://laser.ceb.cam.ac.uk/research/our-software
%
% Please reference:   .......... ??????????? ADD REFERENCE
%
% Notes 
% Co-ordinates
%  XX co-ordinations increase from left to right (Column-co-ord in Matlab)
%  YY co-ordinates increase from top to bottom. (Row co-ord in Matlab)
%  Care is needed in defining orientation and rotation!
% Iterative fit method
%  Starting with Centroid and Orientation guess from regionprops is useful
%  Then optimise by iterative least squares fitting




% 0. SETUP
% Good data for illustrations:
% SUBTILIS:
% fileIn = 'testData_subtilis_ellipsoids.tif'; % SPORE 14 IS A NICE EXAMPLE. Spore 5 for the model data fit plot. 

% SIMULATED OBJECTS:
  fileIn = 'testData_9testSTORMspheres_500nmRad_pixel74nm_PSF145nm.tif'; % TestSTORM data

% MEGATERIUM:
% fileIn = 'testData_megaterium_spheres.tif'; % Spherical megaterium super-resolution       
  


flagModelType = 5;     % 1: Algebraic thin spherical shell 
                       % 2: Monte Carlo thin spherical shell 
                       % 3: Monte Carlo strained spherical shell
                       % 4: Monte Carlo uniformly bright prolate ellipsoid
                       % 5: Monte Carlo ellipsoid, variable equatorial bias

flagOneImageOnly    = 1; % 0: all images. 
                         % 1: One image. 
                         % 2: first N images
                         
flagGetCalledByBatch = 0; % Changes the way this script handles in/out-put
                          % Allows batch script to overwrite the inputs
                          % Causes data to be saved to outFolder, /out
                         
singleImageNumber   = 1; % Index of which single candidate image to fit
flagFirstNimages    = 20; % Or process this many candidates

flagUseTestcardImages = 0; % Use a testcard image of 1 ellipsoid as input
flagSetSegmentManually= 1; % Manually define where to fit model 

flag_WeightBins     = 0;  % Weight pixel data by (1/r) for fitting
flagAntiCollision   = 0;  % Reduce paired candidates to singlets
collisionRadius     = 12; % Anti-collision radius (pixels)

flagShowImages      = 1; % Show image data (in grey)
flagPlotFitOnImages = 1; % Plot model shell outline over image data

flagSaveAllImages   = 0; % Save the images plotted by commands above
flagSaveGoodImages  = 0; % Save the set of accepted images (Not done yet)
flagSaveBadImages   = 0; % Save the set of rejected images (Not done yet)

flagShowCandidates  = 1; % Show raw frame of data with segmentation 
flagShowFittedCenters=1; % Show final fitted centres on raw frame of data
flagShowFitOutlines = 1; % Draw fitted spore outlines over raw data
flagShowTextLabels  = 1; % Show candidate number tags on raw frame of data

flagShowHoffRadHist = 0; % Show histogram of Hough transform rads for debug

flagLimitFitRange   = 1; % Use data in a smaller region of interest:
    lo = 0;              % E.g. for radial positions 0 < r pixels
    hi = 13;             % E.g. -13 < X < 13 and -13 < Y < 13
                         % This is useful to exclude other spores, noise

fL      = 5; % Try 3 - 12 for Florian
fH      = 15;
radSkip = 12; % Try 15 for Florian

flagMedianBGsub     = 0; % Weird background subtraction for research...

fitrad = 20;             % Segment boxes of -fitrad:+fitrad size
                         % Needs to be bigger than "hi" above


% 1. INPUT. 
%    Read in an image
if(flagGetCalledByBatch)
    fileIn = [folder,  batchFileName ];
end

if(flagUseTestcardImages)
 fileIn = ['test.tif'];
 centers = [200,200; 210,210];
 radii   = [8,8];
 metric  = [10,10];
 % flagAntiCollision = 0 % If you want to demonstrate anti-collision filter 
end
if(flagSetSegmentManually)
 centers = [68,72];
 radii   = [4];
 metric  = [10];
end

imDat   = imread(fileIn);
imDatCp  = mean(imDat,3); % Make a grey copy for analysis

if(flagMedianBGsub)
 % imDatCp = 2^16 - 1 - imDatCp;
 fL      = 3; % Try 3 - 12 for Florian
 fH      = 12;
 radSkip = 15; % Try 15 for Florian

 imDatCp = imDatCp - median(imDatCp(:));
 imDatCp(imDatCp < 0) = 0;
end
 
% 2. ANALYSE:
%   (a) Find candidate spores
%   (b) Exclude candidates too near the edge of the camera.

[centers,radii,metric] =imfindcircles(imDatCp,[fL fH],'Sensitivity', 0.90);
 
A = [centers, radii, metric];  % Remove candidates near edge...
  A(A(:,1)<fitrad+8,:) = [];
  A(A(:,1)>size(imDat,2)-(fitrad+8),:) = [];
  A(A(:,2)<fitrad+8,:) = [];
  A(A(:,2)>size(imDat,1)-(fitrad+8),:) = [];
centers = A(:,1:2);
radii   = A(:,3);
metric  = A(:,4);


if(flagAntiCollision)
   lp = 1;
   while lp < length(radii) % For each candidate
      dists = sqrt((centers(:,1)-centers(lp,1)).^2 + (centers(:,2)-centers(lp,2)).^2 );
      dists(lp) = collisionRadius + 100; % Don't exlcude the candidate due to itself
      minDist = min(dists);
      if(minDist<collisionRadius) % Exclude candidate if another is nearby
          centers(lp,:) = [];
          radii(lp) = [];
          metric(lp) = [];
          continue; % Allow list to shorten onto current lp index
      else
      lp = lp + 1;  % Move to next canditate
      end
   end
end

if(flagShowCandidates)
% Plot the camera image with field of candidate spores circled:
figure(1)
  imagesc(imDat);
  colormap(gray)
  truesize; % This is a screensize-dependent bodge for the scatterplot
  hold on
   scatter(centers(:,1),centers(:,2), pi*radii.^2/(9/4),'co','lineWidth',2)
  hold off
  if(flagShowTextLabels)
     a = [1:length(radii)]';
     b = num2str(a); 
     c = cellstr(b);
     dx = 5; 
     dy = -1; % displacement so the text does not overlay the data points
     text(centers(:,1)+dx, centers(:,2)+dy, c, 'color', 'g', 'fontSize', 14);
  end
  if(flagShowHoffRadHist)
  figure(2)
  hist(radii,20)
  xlabel('Segmetation radius, px','fontSize',18)
  ylabel('Number of spores','fontSize',18)
  title('Approx spore size distribution','fontSize',18)
  end
end


% Decide which image data to analyse, and preallocate lists for results
if(flagOneImageOnly ==1)
   startIndex  = singleImageNumber;
   finishIndex = singleImageNumber;
   numberCands = 1;
elseif(flagOneImageOnly ==2)
   startIndex  = 1;
   if(flagFirstNimages < length(radii))
     finishIndex = flagFirstNimages;
     numberCands = flagFirstNimages;
   else
     finishIndex = length(radii);
     numberCands = length(radii);
   end
else
   startIndex  = 1;
   finishIndex = length(radii);
   numberCands = length(radii);
end
    
listFittedRad = zeros(numberCands,1); % Radius, or semi-minor axis length
listFittedAR  = zeros(numberCands,1); % Aspect ratio
listFittedVar = zeros(numberCands,1); % "Sigma" of apparent point spread fn
listFittedRSS = zeros(numberCands,1); % Residual Sum squares, %.
listFittedRow = zeros(numberCands,1); % Row (-ve Y) position in raw data
listFittedCol = zeros(numberCands,1); % Col (+ve X) position in raw data
listFittedInd = startIndex:finishIndex; % Index to match list to image file
listFittedEqu = zeros(numberCands,1); % Equatoriality for biased model
listFittedMax = zeros(numberCands,1); % Maximum brightness in fit
listFittedPsi = zeros(numberCands,1); % Ellipsoid orientation

for lpSpore = startIndex:finishIndex;
    
  cRow = floor(centers(lpSpore,2) );
  cCol = floor(centers(lpSpore,1) );

  imSpore=imDatCp((cRow-fitrad):(cRow+fitrad),(cCol-fitrad):(cCol+fitrad));
  
  % Estimate spore centroid and orientation:
  stats = regionprops(imSpore > mean(imSpore(:)),'Area','Orientation','Centroid');
  areas       = cat(1, stats.Area);
  orientations= cat(1, stats.Orientation);
  centroids   = cat(1, stats.Centroid);
  dat2  = [areas, orientations, centroids];
  dat2  = sortrows(dat2); % Sort by first row (Area) ascending
  orientation = dat2(end,2)*(pi/180); % estimate major axis orientation
  centroid    = dat2(end, 3:4);       % [X, Y] or [COL, ROW] estimate 
  if(abs(centroid(1)-(fitrad+1))>radSkip || abs(centroid(2)-(fitrad+1))>radSkip )
    listFittedRad(lpSpore) = -1;
    listFittedAR(lpSpore)  = -1;
    listFittedVar(lpSpore) = -1;
    continue; % Skip analysis of this spore if 1st guess is way off target    
  end
  
  
  figure(5)
  imagesc(imSpore);
     
  [XX,YY] = meshgrid(1:(2*fitrad+1));
  listX = XX(:);
  listY = YY(:);
  % background = min(imSpore(:));
  background = median(imSpore(imSpore < mean(imSpore(:)))); % 3/15. Good. 
  listI = imSpore(:) - background; % 
  
  % Optionally use only the central part of the image data (less noise):
  if(flagLimitFitRange)
    listI(abs(XX -fitrad - 1) <lo & abs(YY-fitrad - 1)<lo ...
        | abs(XX -fitrad - 1) >hi | abs(YY-fitrad - 1)>hi  ) = [];
    listX(abs(XX -fitrad - 1) <lo & abs(YY-fitrad - 1)<lo ...
        | abs(XX -fitrad - 1) >hi | abs(YY-fitrad - 1)>hi  ) = [];
    listY(abs(XX -fitrad - 1) <lo & abs(YY-fitrad - 1)<lo ...
        | abs(XX -fitrad - 1) >hi  | abs(YY-fitrad - 1)>hi  ) = [];
  end
  
  
  
  % Initial guess of parameters: 
  % [X-centre, Y-centre, radius, sigmaSq, height, ellipticity, azimuth]
  initHt = max(imSpore(:)) - background;
  b0  = [centroid(1),centroid(2),6,9, initHt , 0.2, orientation];
  % radius, sigmaSq, height]
  b0A = [9,9,max(imSpore(:))];
  X = [listX,listY]; % Co-ordinates for the models to use as inputs
  
  if(flagModelType == 1)
      % Exact thin spherical shell model
      mdl = fitnlm(X,listI, @image_DRees, b0(1:5) ); % 
  elseif(flagModelType == 2)
      % Monte Carlo thin spherical shell model
      % mdl = fitnlm(X,listI, @image_sphere_Monte, b0A) ; % Crap!
       beta = fitSphereParams(X, listI, b0(1:6)); % Heuristic least square
  elseif(flagModelType == 3)
      % Monte Carlo strained spherical shell (prolate) model
       beta = fitStrSphParams(X, listI, b0); % Heuristic least squares
  elseif(flagModelType == 4)
      % Monte Carlo uniformly bright prolate ellipsoid model
      beta = fitUnifElParams(X, listI, b0);
  elseif(flagModelType == 5)
      % Monte Carlo prolate ellipsoid model with variable equatorial bias
      Q = 0; % First guess of Equatoriality
      beta = fitBiasElParams(X, listI, [b0, Q]);
  end

  
  if(flagModelType == 1 )   
      mdlC = mdl.Coefficients.Estimate;
   mdlXCen =  mdlC(1); % 
   mdlYCen =  mdlC(2); % 
   mdlRad  =  mdlC(3);
   mdlVar  =  mdlC(4);
   mdlMax  =  mdlC(5);
   mdlEll  =  0; % The sphere has shape factor 1, or (c/a) - 1 = 0;
   mdlPsi  =  0; % And can arbitarily have azimuthal orientation zero. 
   mdlEqu  =  0;
   
   % Evaluate residual sum of squares as a fraction of (sum sq.) image data
   I = image_DRees([mdlXCen,mdlYCen,mdlRad,mdlVar,mdlMax],X);
   sumSq = sum((I - listI).^2); % Difference between exact model and data
   relSumSq = sumSq / sum(listI.^2); % relative sum square error
  elseif(flagModelType == 2)
   mdlXCen =  beta(1);
   mdlYCen =  beta(2);
   mdlRad  =  beta(3);
   mdlVar  =  beta(4);
   mdlMax  =  beta(5);
   mdlEll  =  0;
   mdlPsi  =  0;
   mdlEqu  =  0;
  elseif(flagModelType == 3)
   mdlXCen =  beta(1);
   mdlYCen =  beta(2);
   mdlRad  =  beta(3)
   mdlVar  =  beta(4);
   mdlMax  =  beta(5);
   mdlEll  =  beta(6)
   mdlPsi  =  beta(7);
   mdlEqu  =  -0.5;    % Equatoriality is actually not determined. Is -ve. 
  elseif(flagModelType == 4)
   mdlXCen =  beta(1);
   mdlYCen =  beta(2);
   mdlRad  =  beta(3)
   mdlVar  =  beta(4);
   mdlMax  =  beta(5);
   mdlEll  =  beta(6)
   mdlPsi  =  beta(7); 
   mdlEqu  =  0;
  elseif(flagModelType == 5)
   mdlXCen =  beta(1);
   mdlYCen =  beta(2);
   mdlRad  =  beta(3)
   mdlVar  =  beta(4);
   mdlMax  =  beta(5);
   mdlEll  =  beta(6)
   mdlPsi  =  beta(7);
   mdlEqu  =  beta(8); % Equatoriality
  end
  
  % STORE FITTED RADIUS OF THIS CANDIDATE SPORE IMAGE IN A LIST
  listFittedRad(lpSpore + 1 - startIndex) = mdlRad;
  listFittedAR(lpSpore + 1 - startIndex)  = mdlEll + 1;
  listFittedVar(lpSpore + 1 - startIndex) = mdlVar;
  listFittedRSS(lpSpore + 1 - startIndex) = relSumSq;
  listFittedRow(lpSpore + 1 - startIndex) = cRow + mdlYCen - fitrad - 1;
  listFittedCol(lpSpore + 1 - startIndex) = cCol + mdlXCen - fitrad - 1;
  listFittedEqu(lpSpore + 1 - startIndex) = mdlEqu;
  listFittedMax(lpSpore + 1 - startIndex) = mdlMax;
  listFittedPsi(lpSpore + 1 - startIndex) = mdlPsi;
  
  if(flagShowImages)
     figure(3)
     imagesc(imSpore);
     colormap(gray);
     if(flagPlotFitOnImages)
      hold on
      myPhi = 0:0.01:2*pi; % Draw a circle of correct size on the image...
      if(flagModelType == 1 || flagModelType ==2)    
        myR = mdlYCen(1,1) + mdlRad*cos(myPhi); % SPHERICAL EXACT MODEL
        myC = mdlXCen(1,1) + mdlRad*sin(myPhi);
        plot(myC,myR, 'r','lineWidth', 2);
      elseif(flagModelType == 3 || flagModelType == 4 ...
              || flagModelType == 5) % ELLIPSOID 
        myXp = mdlRad*(mdlEll+1).*cos(myPhi); % Long axis (imagine on X)
        myYp = mdlRad*(1).*sin(myPhi);
        
        % myC, in Col-direction, is new x-coordinate after rotation
        myC  = mdlXCen(1,1) + myXp.*cos(mdlPsi) + myYp.*sin(mdlPsi);
        myR  = mdlYCen(1,1) - myXp.*sin(mdlPsi) + myYp.*cos(mdlPsi);
        plot(myC,myR, 'r','lineWidth', 2);
      end
     hold off
     end  
     if(flagSaveAllImages)
       myIm = getframe(gcf);
       myIm = myIm.cdata;
       % myIm = double(myIm)./255;
       [pathname,filename,extn] = fileparts(fileIn);
       imwrite(myIm, ['C:\Users\user\Documents\Projects\2014_Spores\MATLAB\out\im', filename, int2str(lpSpore) , '.png'])
     end
  end % Finished plotting candidate spore image and model fit
  
end   % Finished fitting model to all candidates

% Remove any fitted values that were set to -1 due to failed analysis:
   listFittedRSS(listFittedRad == -1) = []; % Note "-1" error code in Rad
   listFittedRow(listFittedRad == -1) = [];
   listFittedCol(listFittedRad == -1) = [];
   listFittedInd(listFittedRad == -1) = [];
   listFittedEqu(listFittedRad == -1) = [];
   listFittedMax(listFittedRad == -1) = [];
   listFittedPsi(listFittedRad == -1) = [];
listFittedRad(listFittedRad == -1)   = [];
listFittedAR(listFittedAR == -1)     = [];
listFittedVar(listFittedVar == -1)   = [];


% 3. OUTPUT
%    Plot a histogram of fitted radii (for spheres) or "a" (ellipsoid)
figure(4)
hist(listFittedRad, [0:0.5:15])
   xlim([0 12])
   xlabel('Fitted shell radius, px','fontSize',18)
   ylabel('Number of spores','fontSize',18)

if(flagShowFittedCenters && flagShowCandidates)
  figure(1)
  hold on
    scatter(listFittedCol, listFittedRow, 50, 'rx')
  hold off
end
if(flagShowFitOutlines) 
  % FOR SPHERES:
  if(flagModelType == 1 || flagModelType ==2)
     figure(1)
     hold on
     myPhi = 0:0.005:2*pi;
     for lpP = 1:length(listFittedRad)
       myR = listFittedRow(lpP) + mdlRad*cos(myPhi); % SPHERICAL EXACT MODEL
       myC = listFittedCol(lpP) + mdlRad*sin(myPhi);
       plot(myC,myR, 'r','lineWidth', 2);
     end
     hold off
     axis equal
  end
  % FOR ELLIPSOIDS:
  if(flagModelType == 3 || flagModelType == 4 ...
                        || flagModelType == 5) % ELLIPSOID 
    figure(1)
    hold on
    myPhi = 0:0.01:2*pi;
    for lpP = 1:length(listFittedRad)
       myXp = listFittedRad(lpP).*(listFittedAR(lpP)).*cos(myPhi); % Long
       myYp = listFittedRad(lpP).*sin(myPhi);
        
       % myC, in Col-direction, is new x-coordinate after rotation
       mdlPsi = listFittedPsi(lpP);
       myC  = listFittedCol(lpP) + myXp.*cos(mdlPsi) + myYp.*sin(mdlPsi);
       myR  = listFittedRow(lpP) - myXp.*sin(mdlPsi) + myYp.*cos(mdlPsi);
       plot(myC,myR, 'r','lineWidth', 2);
       axis equal
    end
    hold off
  end
  if(flagShowTextLabels)
     a = [1:length(radii)]';
     b = num2str(a); 
     c = cellstr(b);
     dx = 5; 
     dy = -1; % displacement so the text does not overlay the data points
     hold on
     text(centers(:,1)+dx, centers(:,2)+dy, c, 'color', 'g', 'fontSize', 14);
     hold off
  end
end

% List some (only roughly quality controlled) average results to console
numberCands
estMeanRad = mean(listFittedRad(listFittedRad > 5 & listFittedRad < 14));
estStdRad = std (listFittedRad(listFittedRad > 5 & listFittedRad < 14));
estVar = mean(listFittedVar(listFittedRad > 5 & listFittedRad < 14));
estAsRatio  = mean(abs(listFittedAR(listFittedRad > 5 & listFittedRad < 14)));
mean(listFittedRad(listFittedRad > 5 & listFittedRad < 14))*0.074;
if(flagModelType == 3 || flagModelType == 4) % FOR PROLATE ELLIPSOID of rev 
    estMeanRad = mean(listFittedRad.*(listFittedAR).^(1/3) )
    mean(listFittedRad.*(listFittedAR).^(1/3) )*0.074;
end

if(flagGetCalledByBatch)
    % Save data for each frame of image data during batch analysis
    outFolder = 'C:\Users\user\Documents\Projects\2014_Spores\MATLAB\out\';
    save([outFolder, batchFileName(1:end-4)], ...
           'listFittedRad', 'listFittedRow','listFittedCol', ...
           'listFittedAR', 'listFittedVar', 'listFittedRSS', ...
           'listFittedEqu', 'listFittedMax', ...
           'listFittedInd', 'listFittedPsi', ...
           'estMeanRad', 'estStdRad',...
           'estVar', 'relSumSq' );
end

% THE SHADOW REMAINS CALM!