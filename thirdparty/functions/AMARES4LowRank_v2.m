function DMIResults=AMARES4LowRank_v2(inputfids,Parameters,xaxis,waterloc)
defaultfolder=cd;

%%
% mainfolder = matlab.desktop.editor.getActiveFilename;
% mainfolder = mainfolder(1:size(mainfolder,2)-40);%If the name of the script changes this line has to be modified(-21)
mainfolder = 'E:\mfiles_kmin\MRM_dSPICE\demo\thirdparty\';
addpath(strcat(mainfolder,'\OXSA-master\main'));

%%
AMARESfolder            = strcat(mainfolder,'functions\PriorknowledgeandShowfit');

addpath(AMARESfolder);
DMIpar.samples          = Parameters.NP;
DMIpar.imagingFrequency = Parameters.Freq/(10^6);
DMIpar.timeAxis         = Parameters.time;
DMIpar.dwellTime        = 1/(Parameters.BW);
DMIpar.ppmAxis          = xaxis-4.7;
DMIpar.offset           = -4.7; % This may be modified in the future to include frequency shift based on B0 imperfections
DMIpar.beginTime        = Parameters.TE;

disp(strcat('TE= ',num2str(Parameters.TE*1000),'ms for AMARES'));

% Load prior knowledge file
cd(AMARESfolder)
pk = PK_2H_7T_Simulation; disp('Prior knowledge for Simulation.') % Which prior knowledge will be used is to be decided
% pk = PK_2H_7T_Acetone; disp('Prior knowledge for Acetone.') % Which prior knowledge will be used is to be decided
% pk = PK_2H_7T_LiverDMI_lorentzian;disp('Prior knowledge for Liver DMI.') % Which prior knowledge will be used is to be decided

pkglobal = pk;

for fidnum = 1:numel(inputfids)/DMIpar.samples
    for i = 1:length(pk.initialValues)
        % Adjust prior knowledge values
        DMIpar.offset                 = -(4.7+(4.7-waterloc(fidnum))); % This may be modified in the future to include frequency shift based on B0 imperfections
        pk.initialValues(i).chemShift = (pkglobal.initialValues(i).chemShift + DMIpar.offset);
        pk.bounds(i).chemShift        = (pkglobal.bounds(i).chemShift + DMIpar.offset);
    end

    % Fit with AMARES
    [fitResults, fitStatus, ~, CRBResults] = AMARES.amaresFit(inputfids(:,fidnum), DMIpar, pk, 0);
    [row,col,slice]                        = ind2sub(Parameters.CSIdims,fidnum);
    DMIResults{row,col,slice}              = fitResults;
    DMIResults{row,col,slice}.CRLB         = CRBResults;
    
    cd(defaultfolder)
    fprintf('AMARES fitting progressing... %1.1f%%  \n',fidnum/(numel(inputfids)/DMIpar.samples)*100);
end

disp('AMARES fitting finished.')
end