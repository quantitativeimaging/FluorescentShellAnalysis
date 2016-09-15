% opendata_cylinder_walls
% Demo of cylindrical shell localisation microscopy on sample image data
% 
% Eric Rees 2016 (University of Cambridge)
% CC-BY
% 
% Sample data: Bob Turner, University of Sheffield
% Further work: James Manton, University of Cambridge
%
% OUTLINE
% 
%    This script analyses a pre-defined region of sample image data
%  using cylindrical fluorescent shell localisation microscopy
%    
%    This script generates the cylindrical shell figures for a paper 
%  provisionally called:
%
%   "ELM: super-resolution analysis of wide-field
%   images of fluorescent shell structures"
%
% INSTRUCTIONS
%  1  Ensure this script is in a folder together with the following files:
%      image_cylWall_Monte.m  (forward model for image of cylinder)
%      fitCylWallParams.m     (method to fit parameters of model to data)
%      testData_Bsubtilis168_HADA_cylinders.tif  (sample data)
%  2  Run this script in Matlab
%  
% REFERENCE
%    This method is related to Ellipsoid Localisation Microscopy (ELM), 
%   except the forward model here is for a cylindrical fluorescent shell. 
%    The original ELM paper for spores is published open-access at: 
%   http://www.cell.com/biophysj/fulltext/S0006-3495(15)00993-5  or:
%   http://dx.doi.org/10.1016/j.bpj.2015.09.023 
% 
%   Please cite as appropriate
%
% METHOD
%  1  Read in sample image data
%  2  Select a (pre-defined) region of image data for fitting
%     The region contains a section in the middle of the image of a 
%     cylindrical fluorescent shell, excluding the ends
%  3  Fit cylindrical shell model to image data
%  4  Produce graph and figure output comparing data and fitted model
%  5  Produce super-resolved reconstruction (possibly with overlay).
% 


% 1  INPUT
%    Read in the sample frame of image data

fileIn   = 'testData_Bsubtilis168_HADA_cylinders.tif'; % is RGB uint16 
imDat    = imread(fileIn);
imDatCp  = mean(imDat,3); % Make a grey copy, as double, for analysis

figure(1)
  imagesc(imDatCp);
  colormap(gray);


% 2  ANALYSIS part 1
%    Select a region of image data
%    

% Define a region with corners [xi, yi], pre-chosen by roipoly and copied: 
xi = [263.3909; 280.3820; 271.8864; 253.9513; 263.3909];
yi = [312.9345; 317.0948; 337.0639; 330.4076; 312.9345];

% Create a logical mask for imDatCp from vertices [xi, yi]:
myMask = poly2mask(xi, yi, size(imDatCp,1), size(imDatCp,2));

% Identify x- and y-indexes of each pixel selected in the mask:
[XX,YY] = meshgrid(1:size(imDatCp,2), 1:size(imDatCp,1) );
  listX = XX(myMask);
  listY = YY(myMask);

% For plotting figures using only a selected region of imDatCp:
myBox = [min(listX), min(listY), max(listX)-min(listX), max(listY)-min(listY)];
myROI      = imcrop(imDatCp, myBox );
background = median(myROI(myROI<mean(myROI(:))));
maskBact   = (myROI>mean(myROI(:))); 
mySig      = myROI - background;    % 'signal' may be wanted later
   
listX = listX + 1 - min(listX(:));  % x-indexes in cropped data, myROI
listY = listY + 1 - min(listY(:));

listI   = imDatCp(myMask);
bgFloor = min(listI);               % Floor background keeps listI >=0
listI   = listI - bgFloor;

X       = [listX, listY];           % Put (x,y) coordinates in 1 array
maxdiag = sqrt(max(listX)-min(listX)^2 + (max(listY)-min(listY))^2);

figure(2)
  imagesc(mySig);                        % Show raw image data
  colormap(gray);
  title('Region of image data')
hold on
  plot(xi - min(xi), yi - min(yi), 'g'); % draw region of interest
hold off

figure(3)
  imagesc(maskBact)                      % Show masked region
  title('Mask of approximate bacteria location')


% 3  ANALYSIS part 2
%    Fit a cylindrical fluorescent shell model to the segment

% Obtain a first guess of cylinder position and orientation from the mask
stats = regionprops(maskBact,'Area','Orientation','Centroid');
  areas       = cat(1, stats.Area);
  orientations= cat(1, stats.Orientation);
  centroids   = cat(1, stats.Centroid); % X or col, Y or row

guessX      = centroids(1,1);
guessY      = centroids(1,2);
guessR      = sqrt(areas(1))/3.1;
guessVar    = 5;                   % 'sigma squared' in Gaussian PSF
guessHeight = max(listI(:));
guessPsi    = orientations(1)*(pi/180); % radians (this is far from perfect!)

% Assemble the 1st-guess cylinder shell parameters for fitting
b0 = [guessX, guessY, guessR, guessVar, guessHeight, guessPsi];

% CALL ITERATIVE FITTING METHOD
% Fit cylinder parameters by inverse modelling:
beta = fitCylWallParams(X, listI, b0(1:6)); % Heuristic least square method

% READ OUT SPORE PARAMETERS
mdlXCen =  beta(1);
mdlYCen =  beta(2);
mdlRad  =  beta(3);
mdlVar  =  beta(4);
mdlMax  =  beta(5);
mdlPsi  =  beta(6);


% 4  OUTPUT part 1
%    Produce graph and figure output comparing data and fitted model

% Graph of raw data pixel values versus fitted pixel values 
figure(9)
scatter(ySec, listI, 'bo');
   hold on
     plot(ySec,Ifitted,'rx')
   hold off
   xlabel('Cross section / pixel widths', 'fontSize', 16)
   ylabel('Pixel value', 'fontSize', 16)
set(gca, 'fontSize', 16)
legend('Image data','Model fit')
xlim([-10 10])

% Fitted pixel value shown as an image
imFit = zeros(myBox(4), myBox(3));
  for lp = 1:length(Ifitted)
    imFit(X(lp,2),X(lp,1)) = Ifitted(lp);
  end
figure(11)
  imagesc(imFit); % Note transpose for image orientation!
  colormap(gray)
  caxis([min(myROI(:))-bgFloor, max(myROI(:))-bgFloor])
  title('Fitted model')
hold on
  plot(xi - min(xi), yi - min(yi), 'g'); % draw region of interest
hold off


% 5  OUTPUT part 2
%    Reconstruct super-res. image, based on cylindical shell assumption

%    5(a) Reconstruction INPUT
scaleFc       = 2; % pixel length scale (division) factor  
imSimAll      = zeros(size(imDatCp) * scaleFc ); % Empty image reconstr.n

% Put a simple quality control step here if needed

numberRegions = 1;            % In this script, just deal with one region
numberCyls    = numberRegions; 

lpBatch           = 1;        % Put fitted data in lists, so the reconstruction method could be scaled up
listXi(:,lpBatch) = xi;
listYi(:,lpBatch) = yi;
listRad(lpBatch)  = mdlRad;
listXCen(lpBatch) = mdlXCen + min(listXi(:,lpBatch));
listYCen(lpBatch) = mdlYCen + min(listYi(:,lpBatch));
listVar(lpBatch)  = mdlVar;
listMax(lpBatch)  = mdlMax;
listPsi(lpBatch)  = mdlPsi;
listDiags(lpBatch)= maxdiag;
listInd(lpBatch)  = lpBatch;

imOverlaid = imresize(imDatCp, scaleFc, 'nearest'); % No-interpolation scale-up of image data

%    5(b) Reconstruction Process
for lpCy = 1:numberCyls
 indexParams = listInd(lpCy); % Some may have been deleted by qual. contr.

 sRow = floor(listYCen(listInd == indexParams)); % Bodged method to read data for one region
 sCol = floor(listXCen(listInd == indexParams));
 
 params = [listXCen(listInd == indexParams), ... 
           listYCen(listInd == indexParams), ...
           listRad(listInd == indexParams), ...
           0.15, ... % Variance = applied PSF-like blur, 0.15 ~~ 0.4^2
           listMax(listInd == indexParams), ...
           listPsi(listInd == indexParams), ...
           ];
    
 paramsScaled = [params(1:4)*scaleFc, params(5:end)];
 
 maxdiag = listDiags( listInd==indexParams );
 fitrad  = ceil(maxdiag)*1.0;
 
 % Make a rectangular meshgrid of pixel co-ordinates in the reconstruction
 [XX,YY] = meshgrid( (sCol*scaleFc-fitrad):(sCol*scaleFc+fitrad),...
                     (sRow*scaleFc-fitrad):(sRow*scaleFc+fitrad) );
                
 % Hence produce a list of (x,y) co-ordinates, [listX, listY]
 % To only overlay reconstruction inside the fitting polygon region
 flagLimitReconToPolygons = 1;      
 if(flagLimitReconToPolygons  == 1)
   [XXX,YYY] = meshgrid(1:size(imDatCp,2)*scaleFc, 1:size(imDatCp,1)*scaleFc);
   
   xip = listXi(:,lpCy)*scaleFc;
   yip = listYi(:,lpCy)*scaleFc;
   sz1 = size(imDatCp,2)*scaleFc;
   sz2 = size(imDatCp,1)*scaleFc;
   
   bMask = poly2mask(xip, yip, size(imDatCp,1)*scaleFc, size(imDatCp,2)*scaleFc);

   listX = XXX(bMask);
   listY = YYY(bMask);
 else
   listX = XX(:);
   listY = YY(:);
 end
 X = [listX,listY];
 
 % Simulate pixel values using the cylindrical image model
 I = image_cylWall_Monte(paramsScaled, X);
 
 % Add pixel values to reconstuction
 for lp = 1:length(I)
   imSimAll(X(lp,2),X(lp,1)) = imSimAll(X(lp,2),X(lp,1)) + I(lp);
 end
 
 % leave the following commented to not re-add d.c. background:
 % imSimAll = imSimAll + bgFloor;


for lp = 1:length(I)
 imOverlaid(X(lp,2)-1,X(lp,1)) = I(lp);
end

 
end

% Crop a local region of imOverlaid for more compact comparison
try
  rect = [min(xip)-8, min(yip)-8,max(xip)-min(xip)+16 max(yip)-min(yip)+16];
  imSegment = imcrop(imOverlaid, rect);
catch ME
  warn('Box for cropping sub-recon goes outside bounds of recon');
end
  
%    5(c) Reconstruction output
figure(12)
imagesc(imSimAll)
colormap(gray)
axis equal
%  hold on
%    scatter(listXCen*scaleFc,listYCen*scaleFc,'r+')
%  hold off

figure(13)
  imagesc(imOverlaid)
  colormap(gray)
axis equal

figure(14)
  imagesc(imSegment)
  colormap(gray)
axis equal

%    For image-saving convenience:
flagSaveImages = 0;
if(flagSaveImages)
     imwrite(imSimAll./max(imSimAll(:)), 'simIm.png')

     imwrite(imOverlaid./max(imSimAll(:)), 'simOverlay.png')
end

