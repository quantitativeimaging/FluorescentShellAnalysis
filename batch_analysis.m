% Script to run spore image analysis on many image stacks
% 
% Notes
%   Remember to check this script calls the correct version of the fitting
%  script! E.g. the current version is v5:  spores_analysis_v5_ellipsoids;


folders = {'C:\Users\user\Documents\Projects\2014_Spores\2015_subtilis2\'};

% {'C:\Users\user\Documents\Projects\2014_Spores\2015_April_Subtilis GFP spores\Subtilis GFP spores\', ...
%     'C:\Users\user\Documents\Projects\2014_Spores\2015_March\GFP Images for Eric\Inner Coat\SleL\', ...
%     'C:\Users\user\Documents\Projects\2014_Spores\2015_March\GFP Images for Eric\Inner Coat\3035\', ...
%     'C:\Users\user\Documents\Projects\2014_Spores\2015_March\GFP Images for Eric\Inner Coat\4051\',...
%     'C:\Users\user\Documents\Projects\2014_Spores\2015_March\GFP Images for Eric\Outer Coat\CotE\',...
%     'C:\Users\user\Documents\Projects\2014_Spores\2015_March\GFP Images for Eric\Outer Coat\737\',...
%     'C:\Users\user\Documents\Projects\2014_Spores\2015_March\GFP Images for Eric\Exosporium\pHT313CotX28GFP in PV361\'...
%     'C:\Users\user\Documents\Projects\2014_Spores\2015_March\GFP Images for Eric\Exosporium\CotX29_FL3_GFP\'...
%     'C:\Users\user\Documents\Projects\2014_Spores\18Dec2014\Super resolution Images\CotX28 QMB\', ...
%     'C:\Users\user\Documents\Projects\2014_Spores\2015_March\GFP Images for Eric\Exosporium\CotW\'};
    

for lpFolder = 1:size(folders,2)
folder = folders{lpFolder};

 listOfFiles = dir(folder);

 numberOfBatchFiles = length(listOfFiles);

 for lpBatch = 1:numberOfBatchFiles

    batchFileName = listOfFiles(lpBatch).name;
    if(batchFileName(1)=='.')
       continue % Skip hidden files
    end
   
    spores_analysis_v5_ellipsoids;
 end

end