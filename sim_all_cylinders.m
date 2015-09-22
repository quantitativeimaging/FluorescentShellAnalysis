% sim_all_cylinders
% 
% Generates a reconstructed image of several cylindrical fluorescent shells
% Based on parameters that are obtained by the cylindrical shell analysis
% 
% Notes
%
%
% To use this script:
%  1. First, run the analysis script, batch_cylWall.m, on sample data 
%     such as testData_Bsubtilis168_HADA_cylinders.tif
%  2. Then run this script, which will read parameters from the base 
%     workspace.


% 0. SETUP
flagSaveImage           = 0;
flagCaptureListedParams = 1;
% flagQualControl       = 1; % Apply some quality control step?
%   cutVar              = 30;
% modelType             = 1;   % Select image model. I've 1 cylinder model.

% Parameters for reconstruction
fitrad                  = 24;  % 
% nPoints            = 10000;  % Cylinder wall model lacks this input
scaleFc                 = 2;   % Pixel scale factor for reconstruction


% 1. INPUT
imSim      = zeros(size(imDatCp) * scaleFc ); % Empty image reconstr.n

% Put a simple quality control step here if needed

numberCyls = numberRegions;    % In case quality control voids some regions

% 2. RECONSTRUCTION
for lpCy = 1:numberCyls
 indexParams = listInd(lpCy); % Some may have been deleted by qual. contr.

 sRow = floor(listYCen(listInd == indexParams)); % Really?
 sCol = floor(listXCen(listInd == indexParams));
 
 params = [listXCen(listInd == indexParams), ... 
           listYCen(listInd == indexParams), ...
           listRad(listInd == indexParams), ...
           0.15, ... % Variance = applied PSF-like blur
           listMax(listInd == indexParams), ...
           listPsi(listInd == indexParams), ...
           ];
    
 paramsScaled = [params(1:4)*scaleFc, params(5:end)];
 
 maxdiag = listDiags( listInd==indexParams );
 fitrad  = ceil(maxdiag)*1.0;
 
 % Generate the local meshgrid of pixel co-ordinates in the reconstruction
 [XX,YY] = meshgrid( (sCol*scaleFc-fitrad):(sCol*scaleFc+fitrad),...
                    (sRow*scaleFc-fitrad):(sRow*scaleFc+fitrad) );
 % Hence produce a list of (x,y) co-ordinates:
 listX = XX(:);
 listY = YY(:);
 X = [listX,listY];
 
 % Simulate pixel values from image model
 I = image_cylWall_Monte(paramsScaled, X);
 
 % Add pixel values to reconstuction
 for lp = 1:length(I)
   imSim(X(lp,2),X(lp,1)) = imSim(X(lp,2),X(lp,1)) + I(lp);
 end

end
    
% 3. OUTPUT
figure(10)
imagesc(imSim)
colormap(gray)
axis equal
%  hold on
%    scatter(listXCen*scaleFc,listYCen*scaleFc,'r+')
%  hold off

%    For image-saving convenience:
if(flagSaveImage)
     imwrite(imSim'./max(imSim(:)), 'simIm.png')
end