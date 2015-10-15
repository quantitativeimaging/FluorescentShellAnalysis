% sim_all
%
%  Simulate the images of all the spores based on their assessed parameters
%  This produces a reconstructed image rather similar in principle to a 
% fluorescence density reconstruction in localisation microscopy.
%  The reconstructed image will only be a realistic image under the 
% assumption that the model correctly describes the spore shells and also 
% that the analysis correctly captures the shell parameters
%

flagSaveImage           = 0;
flagCaptureListedParams = 1;

flagQualControl = 1;
  cutVar = 30;
  
modelType = 1; 

% Size in camera pixels
fitrad = 24;
nPoints = 10000;

scaleFactor = 2;

imSim = zeros(fliplr(size(imDatCp)) * scaleFactor );


% Quality control: filter out Var >> expected
if(flagQualControl)
  listFittedInd = listFittedInd(listFittedVar < cutVar);
  listFittedCol = listFittedCol(listFittedVar < cutVar);
  listFittedRow = listFittedRow(listFittedVar < cutVar);
  listFittedRad = listFittedRad(listFittedVar < cutVar);
  listFittedMax = listFittedMax(listFittedVar < cutVar);
  listFittedAR  = listFittedAR(listFittedVar < cutVar);
  listFittedPsi = listFittedPsi(listFittedVar < cutVar);
  listFittedEqu = listFittedEqu(listFittedVar < cutVar);
  listFittedVar = listFittedVar(listFittedVar < cutVar);
end
  

for lpSpore = 1:length(listFittedInd)
 indexParams = listFittedInd(lpSpore) % First on list (!) or, say 45

 sRow = floor(listFittedRow(listFittedInd == indexParams));
 sCol = floor(listFittedCol(listFittedInd == indexParams));
 
 params = [listFittedCol(listFittedInd == indexParams), ...
           listFittedRow(listFittedInd == indexParams), ...
           listFittedRad(listFittedInd == indexParams), ...
           1, ... % Variance = applied PSF-like blur  
           listFittedMax(listFittedInd == indexParams), ...
           listFittedAR(listFittedInd == indexParams) - 1, ...
           listFittedPsi(listFittedInd == indexParams), ...
           listFittedEqu(listFittedInd == indexParams), ...
           ];

  params2 = [params(1:4)*scaleFactor, params(5:8)];
  % imSim = zeros(2*fitrad*scaleFactor + 1)
  %[XX,YY] = meshgrid(1:(2*fitrad*scaleFactor+1));
  [XX,YY] = meshgrid( (sCol*scaleFactor-fitrad):(sCol*scaleFactor+fitrad),...
                    (sRow*scaleFactor-fitrad):(sRow*scaleFactor+fitrad) );
  listX = XX(:);
  listY = YY(:);
  X = [listX,listY];
  
  if(modelType == 1)
    I = image_DRees(params2, X);
  elseif(modelType == 5)
    I = image_biasEl_Monte([params2, nPoints], X);
  end

  for lp = 1:length(I)
    imSim(X(lp,1),X(lp,2)) = imSim(X(lp,1),X(lp,2)) + I(lp);
  end
 
end

figure(9)
imagesc(imSim')
 
colormap(gray)
axis equal
 
%
if(flagSaveImage)
     imwrite(imSim'./max(imSim(:)), 'simIm.png')
end