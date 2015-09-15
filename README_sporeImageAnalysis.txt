% SPORE IMAGE ANALYSIS SCRIPTS.
%
% DETERMINES Fluorescent shell radius and/or other parameters
% by fitting models to image data.
% Uses a segmentation, anti-collision, and fitting method.
% A quality control criteria may be needed to select accurate fits. 
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
% Please reference:   
%  Ellipsoid localisation microscopy infers the size and order of protein layers in Bacillus spore coats
%  Biophys Journal
%
%  J. Manetsberger,1 J.D. Manton,1 M.J Erdelyi,2 H. Lin,1 H.D. Rees, G. Christie,1 E.J. Rees1* 
%  1 Department of Chemical Engineering and Biotechnology, University of Cambridge, CB2 3RA, UK.
%  2 University of Szeged, Department of Optics and Quantum Electronics, Szeged, Dom ter 9, Hungary. 
%  (Submitted - contact ejr36@cam.ac.uk)


Brief summary: 

Capabilities:
1. Five models: (1) exact sphere, Monte Carlo: (2) sphere; (3) stretched sphere (prolate ellipsoid of revolution); (4) uniform ellipsoid; (5) biased ellipsoid (uniform model with added parameter for polar or equatorial bias of the shell thickness)
2. Can output tabulated radii and other parameters in batch mode
3. Can output residual errors to find best model
4. Can generate videos of fitting process


Outputs:
listFittedRad % Sphere radius, or semiminor axis of ellipsoid
listFittedAR  % Aspect ratio (b/a) of the ellipsoid. 1 for sphere models. 
listFittedVar % Fitted variance "sigma squared" of PSF 
listFittedRSS % Residual sum of squares (misfit). Often better than 3%!
listFittedRow % Fitted Row (negative Y) position
listFittedCol % Fitted Col (positive X) position
listFittedInd % Index of spore for reference to raw image data


Notes:
SEGMENTATION
 Circular Hough Transform finds spheroidal spore images. 
 It also finds ellipsoidal images, but often finds both poles separately
 Hence anti-collision filtering to remove duplicated segments

FITTING
 For equation, use standard library function
 For Monte Carlo models, a custom iterative least squares method is used.

QUALITY CONTROL CRITERION.
A "plausible variance" criterion for quality control may be useful.
Selecting only candidates with listFittedVar < 16 (about 300 nm PSF standard deviation in sample data) corresponds to fluorescent shells with a small image blur.
This often accepts all individual spores, while rejecting clumps of spores and fragments of fluorescent matter.

GETTING STARTED: 
1. Start Matlab
2. Navigate to this folder, and run the script "" 
3. It will automatically run a quick analysis on the sample data included in this folder
--- in order to fit models to the full dataset, you will need to edit parameters in the code
--- Try changing to:

flagModelType = 1;       %  Fit the spherical model instead of the ellipsoid (much faster)
flagOneImageOnly    = 0; % 0: Process all images. 


