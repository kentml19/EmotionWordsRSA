% Recipe_fMRI_searchlight
%
% Cai Wingfield 11-2009, 2-2010, 3-2010, 8-2010
%__________________________________________________________________________
% Copyright (C) 2010 Medical Research Council

%%%%%%%%%%%%%%%%%
%% Searchlight %%
%%%%%%%%%%%%%%%%%
clear;clc;
toolboxRoot = 'C:/Users/kentml19/Desktop/Research/Postdoc/1_fMRI/rsatoolbox-develop'; addpath(genpath(toolboxRoot)); % Catch sight of the toolbox code
addpath(genpath(toolboxRoot));
mkdir('EmoWordsRSA');

%Generate userOptions
%userOptions = EmoWordsUserOptions();
userOptions = defineUserOptions();
userOptions.rootPath = toolboxRoot;
userOptions.analysisName = 'EmoWordsRSA';

% Generate a simulationOptions structure.
searchlightOptions.monitor = false;
searchlightOptions.fisher = true;
searchlightOptions.nSessions = 1;
searchlightOptions.nConditions = 10;
models = rsa.constructModelRDMs(EmoWords_modelRDMs, userOptions);

%binaryMasks_nS = rsa.fmri.fMRIMaskPreparation(userOptions);

nCond = userOptions.conditionLabels;
Nsubjects = 15;
%% Grabbing Betas for RSA from NeuroElf (Note run separately for different runs)

y = xff; %Assigns loaded objects in neuroelf to variable
glm = y(5); %Number in parenthesis is object in list (provided by variable 'y' or function 'xff')
vmp = y(7); %See above VMP = statsmap
vmr = y(4); %see above, VMR = anatomical map
mask = vmp.Map.VMPData;
mask(:,:,:) = 1;

%% Searchlight RDMs
for subI = 1:Nsubjects
   thisSubject = userOptions.subjectNames{subI};
   for session = 1
  betas = glm.GLMData.Subject(subI).BetaMaps(:,:,:,:);
  betas(:,:,:,11) = []; %remove the constant regressor
   for condition = 1:10 %10 regressors
       betasW = betas(:,:,:,condition);
       brainMatrix = betasW;
        brainVector = reshape(brainMatrix, 1, []);
				subjectMatrix(:, condition, session) = brainVector; % (voxel, condition, session)
                
                clear brainMatrix brainVector
       end
       subName = userOptions.subjectNames{subI};
    end

                fullBrainVols = subjectMatrix; clear subjectMatrix
                
    subject = ['subject',num2str(subI)];
    maskName = 'mask';
    userOptions.searchlightRadius = 9;
    fprintf(['computing correlation maps for subject %d \n'],subI)
    [rs, ps, ns, searchlightRDMs.(subject)] = rsa.fmri.searchlightMapping_fMRI(fullBrainVols, models, mask, userOptions, searchlightOptions);
    rsa.util.gotoDir(userOptions.rootPath, 'Maps3');
    save(['rs_',subName,'_2','.mat'],'rs');
    clear rs searchlightRDMs;
    clear subName
    cd(toolboxRoot);
    
end

%% load the previously computed rMaps and concatenate across subjects
% prepare the rMaps:
Nsubjects = 15;
for subI = 1:Nsubjects
    subName = userOptions.subjectNames{subI};
    %load([userOptions.rootPath,filesep,'Maps',filesep,'rs_',subName,'.mat']);
    load('EmoWordsRSA_fMRISearchlight_Maps.mat');
%    rMaps = rs;
    rMaps = rMaps_nS.allSeparate.(subName).mask;
    fprintf(['loading the correlation maps for subject %d \n'],subI);
end
% concatenate across subjects
for modelI = 1:numel(models)
    for subI = 1:Nsubjects
        thisRs = rMaps(subI);
        %mask = x;
        %mask = binaryMasks_nS.(subName).mask;
        thisModelSims(:,:,:,subI) = thisRs(:,:,:,modelI);
    end
    % obtain a pMaps from applying a 1-sided signrank test and also t-test to
    % the model similarities:
    for x=1:size(thisModelSims,1)
        for y=1:size(thisModelSims,2)
            for z=1:size(thisModelSims,3)
                if mask(x,y,z) == 1
                    [h p1(x,y,z)] = ttest(squeeze(thisModelSims(x,y,z,:)),0,0.05,'right');
                    [p2(x,y,z)] = rsa.stat.signrank_onesided(squeeze(thisModelSims(x,y,z,:)));
                else
                    p1(x,y,z) = NaN;
                    p2(x,y,z) = NaN;
                end
            end
        end
        disp(x);
    end
    % apply FDR correction
    pThrsh_t  = rsa.stat.FDRthreshold(p1,0.05,mask);
    pThrsh_sr = rsa.stat.FDRthreshold(p2,0.05,mask);
    p_bnf = 0.05/sum(mask(:));
    % mark the suprathreshold voxels in yellow
    supraThreshMarked_t = zeros(size(p1));
    supraThreshMarked_t(p1 <= pThrsh_t) = 1;
    supraThreshMarked_sr = zeros(size(p2));
    supraThreshMarked_sr(p2 <= pThrsh_sr) = 1;
    
    % display the location where the effect was inserted (in green):
    brainVol = rsa.fmri.addRoiToVol(rsa.util.map2vol(vmr.VMRData),rsa.util.mask2roi(mask),[1 0 0],2);
    brainVol_effectLoc = rsa.fmri.addBinaryMapToVol(brainVol,mask.*mask,[0 1 0]);
    rsa.fig.showVol(brainVol_effectLoc,'simulated effect [green]',2);
    rsa.fig.handleCurrentFigure([toolboxRoot,filesep,'EmoWordsRSA',filesep,'results_DEMO4_simulatedEffectRegion'],userOptions);
    
    % display the FDR-thresholded maps on a sample anatomy (signed rank test) :
    brainVol = rsa.fmri.addRoiToVol(rsa.util.map2vol(vmr.VMRData),rsa.util.mask2roi(mask),[1 0 0],2);
    brainVol_sr = rsa.fmri.addBinaryMapToVol(brainVol,supraThreshMarked_sr.*mask,[1 1 0]);
    rsa.fig.showVol(brainVol_sr,'signrank, E(FDR) < .05',3)
    rsa.fig.handleCurrentFigure([toolboxRoot,filesep,'EmoWordsRSA',filesep,'results_DEMO4_signRank'],userOptions);
    
    % display the FDR-thresholded maps on a sample anatomy (t-test) :
    brainVol = rsa.fmri.addRoiToVol(rsa.util.map2vol(vmr.VMRData),rsa.util.mask2roi(mask),[1 0 0],2);
    brainVol_t = rsa.fmri.addBinaryMapToVol(brainVol,supraThreshMarked_t.*mask,[1 1 0]);
    rsa.fig.showVol(brainVol_t,'t-test, E(FDR) < .05',4)
    rsa.fig.handleCurrentFigure([toolboxRoot,filesep,'EmoWordsRSA',filesep,'results_DEMO2_tTest'],userOptions);
end

cd(toolboxRoot);
