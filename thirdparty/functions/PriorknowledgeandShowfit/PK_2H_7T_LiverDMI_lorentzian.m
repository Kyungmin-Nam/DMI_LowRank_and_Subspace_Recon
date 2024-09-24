function outStruct = PK_2H_7T_LiverDMI_lorentzian()
%% .M file to assemble the bounds, priorKnowledge and initialValues structs for the matlab implementation of AMARES
%prior knowledge for in vivo liver 7T compounds
%created by Ayhan Gursan

% No gaussian part included
%%
%Each of B, PK and IV is a 1xN struct, where N is the number of peaks. Note
%multiplets are counted as one peak.
%The fields are as follows:
%bounds           initialValues          priorKnowledge

%peakName         peakName               peakName
%chemShift        chemShift              multiplet
%linewidth        linewidth              chemShiftDelta
%amplitude        amplitude              amplitudeRatio
%phase            phase                  G_linewidth
%chemShiftDelta                          G_amplitude
%amplitudeRatio                          G_phase
%                                        G_chemShiftDelta
%                                        refPeak


%% Bounds
fields.Bounds = {
    'peakName',                                 'chemShift',     'linewidth',   'amplitude',    'phase',     'chemShiftDelta',   'amplitudeRatio'};
values.boundsCellArray = {...
    '1.3Lipid',                               [1.10,1.50],       [10,50],        [0,inf],     [0,360],      [],                [];
    '2.25GLX',                                [2.05,2.45],       [10,50],        [0,inf],     [0,360],      [],                [];
    '3.8Glucose',                             [3.70,3.90],       [10,50],        [0,inf],     [0,360],      [],                [];
    '4.7water',                               [4.65,4.75],       [10,50],        [0,inf],     [0,360],      [],                [];
    };

%% initialValues
fields.IV = {
    'peakName',                                   'chemShift',     'linewidth',   'amplitude',    'phase'};
values.IVCellArray = {...
    '1.3Lipid',                                         1.30,               25,        3.00,               0;
    '2.25GLX',                                          2.25,               25,        3.00,               0;
    '3.8Glucose',                                       3.80,               25,        3.00,               0;
    '4.7water',                                         4.70,               25,        3.00,               0;
    };

%%
fields.PK = {
    'peakName',                                 'multiplet',     'chemShiftDelta',   'amplitudeRatio',    'G_linewidth',   'G_amplitude',    'G_phase'     'G_chemShiftDelta',   'refPeak'};
values.PKCellArray = {...
    '1.3Lipid',                                       [],             [],               [],                  [],            [],               [1],           [],                    0;
    '2.25GLX',                                        [],             [],               [],                  [],            [],               [1],           [],                    0;
    '3.8Glucose',                                     [],             [],               [],                  [],            [],               [1],           [],                    0;
    '4.7water',                                       [],             [],               [],                  [],            [],               [1],           [],                    0;
    };

%% Pass to the function which assembles the constraints into structs and saves them
outStruct = AMARES.priorKnowledge.preparePriorKnowledge(fields,values);
outStruct.svnVersion = '$Rev: 8177 $';
outStruct.svnHeader = '$Header: https://cardiosvn.fmrib.ox.ac.uk/repos/crodgers/FromJalapeno/MATLAB/RodgersSpectroToolsV2/main/+AMARES/+priorKnowledge/PK_H_7T_invivoLiver.m 8177 2016-08-12 km $';
