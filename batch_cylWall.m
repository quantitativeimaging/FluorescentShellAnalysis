% batch_cylWall
%  A script to process several user-selected regions of cylindrical shells
%
% Eric Rees 2015
% License: GNU LGPL version 3+

% INSTRUCTIONS:
% 
%    Move the scripts (batch_cylWall, wall_analysis_v1, image_cylWall_Monte 
%   and fitCylWallMonte) to the same folder.
%
% 1. In wall_analysis_v1.m make sure flagGetCalled = 1; is set
% 2. In this script, batch_cylWall.m, name the input file, including .tif
% 3. In this script, set a suitable number of regions to process
% 4. Run this script
% 5. When Figure 1 appears, use the cursor to select regions of interest
%     Choose straight, cylindrical segments of specimens
%     Select the feature and some surrounding dark space
%     Try to select a region with greater length along the specimen than
%     its width, otherwise the auto-initial guess of orientation may be bad
% 6. Wait for the inverse modelling to process (may take 30 s per segment)
% 7. Multiply listRads by pixel width to get physical cylinder radius
% 8. A quality control step may be needed to reject misfitted regions. 
%     


% 1. INPUT
suppliedFileIn = ['testData_Bsubtilis168_HADA_cylinders']; % 

numberRegions = 1; 


% 2. ESTABLISH WHICH REGIONS TO PROCESS
imDat   = imread(suppliedFileIn);
imDatCp  = mean(imDat,3); % Make a grey copy for analysis

figure(1)
  imagesc(imDatCp);
  colormap(gray);
 
listMasks = false([size(imDatCp),numberRegions]); %
listXi    = zeros(5, numberRegions);
listYi    = zeros(5, numberRegions);
listRad   = zeros(numberRegions, 1);
listXCen   = zeros(numberRegions, 1);
listYCen   = zeros(numberRegions, 1);
listVar   = zeros(numberRegions, 1);
listMax   = zeros(numberRegions, 1);
listPsi   = zeros(numberRegions, 1);

dx = 4; 
dy = 1; % displacement so the text does not overlay the data points

for lpUser = 1:numberRegions
    
    [aMask, xi, yi] = roipoly;
    listMasks(:,:,lpUser) = aMask;
    if(length(xi)>=5) % If a quadrilateral or higher is chosen
      listXi(:,lpUser)     = xi(1:5); % record 4 vertices. 
      listYi(:,lpUser)     = yi(1:5);
    end
    figure(1)
    hold on
     c = int2str(lpUser);
     text(max(xi)+dx, max(yi)+dy, c, 'color', 'g','fontSize',14);
     plot(xi, yi, 'g')
    hold off
    
end
% 3. CALL WALL ANALYSIS AND SAVE FITTED PARAMETERS

for lpBatch = 1:numberRegions
 
    suppliedMask = listMasks(:,:,lpBatch);
    
    xi = listXi(:,lpBatch);
    yi = listYi(:,lpBatch);
    
    wall_analysis_v1;
    
    listRad(lpBatch)  = mdlRad;
    listXCen(lpBatch) = mdlXCen;
    listYCen(lpBatch) = mdlYCen;
    listVar(lpBatch)  = mdlVar;
    listMax(lpBatch)  = mdlMax;
    listPsi(lpBatch)  = mdlPsi;
  
end