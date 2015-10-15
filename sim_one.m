% sim_one
%
% Simulate the image of a spore based on a model and defined parameters
% This can be used to generate something approximately like a super-
% resolution image, by setting the diffraction-limited PSF to a small value

flagSaveImage = 1;

flagModelType = 6;

flagCaptureListedParams = 0;

% Size in camera pixels
fitrad = 20;
nPoints = 30000;


scaleFactor = 1;

% Parameters in pixels. 6.48, 7.84
%          xcen,     ycen, rad or a,  var, height, b/a -1, Psi, equ
% params = [fitrad+1,fitrad+1, 5.10,   1, 991,    0.59,   0, -0.12 ];%SleL
% params = [fitrad+1,fitrad+1, 5.44,   1, 991,    0.70,   0, -0.06 ];% CotZ
%params = [fitrad+1,fitrad+1, 5.56,   1, 991,    0.57,   0, -0.08 ];% CotG
% params = [fitrad+1,fitrad+1, 5.82,   1, 900,    0.50,   0, +0.25 ];%test
 
% example meg sim for one spore: the one shown in the paper
% parameters via load: example_data_model_sphere in paper /  2015/ Matlab
 params = [fitrad+1,fitrad+1, 10.0,   5^2, 900,    0,   0, +0, 7 ]; 

% Example for comparing finite wall thickness results
 
%
% B. meg
% params = [fitrad+1,fitrad+1, 7.51,   1, 991,    0,   0, -0.0 ];%SleL
% params = [fitrad+1,fitrad+1, 7.10,   1, 991,    0,   0, -0.0 ];%3035
% params = [fitrad+1,fitrad+1, 7.67,   1, 991,    0,   0, -0.0 ];%CotE
% params = [fitrad+1,fitrad+1, 8.16,   1, 991,    0,   0, -0.0 ];%CotW

if(flagCaptureListedParams)
   indexParams = listFittedInd(1); % First on list (!) or, say 45
    
    params = [fitrad+1,fitrad+1, ...
              listFittedRad(listFittedInd == indexParams), ...
              1, ... % Variance = applied PSF-like blur  
              listFittedMax(listFittedInd == indexParams), ...
              listFittedAR(listFittedInd == indexParams) - 1, ...
              listFittedPsi(listFittedInd == indexParams), ...
              listFittedEqu(listFittedInd == indexParams), ...
              ];
end

params2 = [params(1:4)*scaleFactor, params(5:8)];

imSim = zeros(2*fitrad*scaleFactor + 1);

[XX,YY] = meshgrid(1:(2*fitrad*scaleFactor+1));
  listX = XX(:);
  listY = YY(:);
X = [listX,listY];

if(flagModelType == 1)
  I = image_DRees([params2, nPoints], X);
elseif(flagModelType == 5)
  I = image_biasEl_Monte([params2, nPoints], X);
elseif(flagModelType ==6)
  I = image_ERees(params, X);
end

for lp = 1:length(I)
    imSim(X(lp,1),X(lp,2)) = I(lp);
end
 
 figure(9)
 imagesc(imSim')
  axis equal
 colormap(gray)
 
 %
 if(flagSaveImage)
     imwrite(imSim'./max(imSim(:)), 'simIm.png')
 end
 %