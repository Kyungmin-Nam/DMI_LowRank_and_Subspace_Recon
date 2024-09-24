function outStruct = PK_2H_7T_Acetone()
%% .M file to assemble the bounds, priorKnowledge and initialValues structs for the matlab implementation of AMARES
%prior knowledge for 7T deuteriated compounds in Acetone Phantom
%created by Ayhan Gursan

% Upper bound for linewidths are set to 10= minimum 32ms of T2star
% Less weighting on lorentzian part

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
    'peakName',                                 'chemShift',    'linewidth',   'amplitude',    'phase',     'chemShiftDelta',   'amplitudeRatio', 'sigma'};
values.boundsCellArray = {...
    '2.4Acetone',                                [2.30,2.50],       [5,40],       [0,inf],     [0,360],      [],                [],       [];
    '4.7water',                               [4.60,4.75],       [5,40],       [0,inf],     [0,360],      [],                [],       [];
    };

%% initialValues
fields.IV = {
    'peakName',                                   'chemShift',    'linewidth', 'amplitude',    'phase',  'sigma'};
values.IVCellArray = {...
    '2.4Acetone',                                          2.40,      20,              3.00,               0,   0;
    '4.7water',                                         4.70,      20,                3.00,               0,   0;
    };

%%
fields.PK = {
    'peakName',                                 'multiplet',     'chemShiftDelta',   'amplitudeRatio',    'G_linewidth',   'G_amplitude',    'G_phase'     'G_chemShiftDelta',   'refPeak'};
values.PKCellArray = {...
    '2.4Acetone',                                        [],             [],               [],                  [],            [],               [1],           [],                    0;
    '4.7water',                                       [],             [],               [],                  [],            [],               [1],           [],                    0;
    };

%% Pass to the function which assembles the constraints into structs and saves them
outStruct = AMARES.priorKnowledge.preparePriorKnowledge(fields,values);
outStruct.svnVersion = '$Rev: 8177 $';
outStruct.svnHeader = '$Header: https://cardiosvn.fmrib.ox.ac.uk/repos/crodgers/FromJalapeno/MATLAB/RodgersSpectroToolsV2/main/+AMARES/+priorKnowledge/PK_H_7T_invivoLiver.m 8177 2016-08-12 km $';
