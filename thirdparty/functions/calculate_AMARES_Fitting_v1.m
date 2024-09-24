function SimPhantom = calculate_AMARES_Fitting_v1(SimParam,ksp)

%--------------------------------------------------------------------------
% image space (x-y-z) <=> k-space (kx-ky-kz)
%--------------------------------------------------------------------------
img_t = ksp;
for dim = 1:3
    img_t = sqrt(size(img_t,dim)) * ...
        fftshift(ifft(ifftshift(img_t, dim), [], dim), dim);
end

% swithch the dimension order of data
SimPhantom.img_t = permute(squeeze(double(img_t)), [4 1 2 3]);

%% Set parameters
SimPhantom.Param.dims               = size(permute(SimPhantom.img_t, [2 3 4 1]));
SimPhantom.Param.Freq               = SimParam.Freq;
SimPhantom.Param.BW                 = SimParam.BW;
SimPhantom.Param.time               = SimParam.time;
SimPhantom.Param.TE                 = 0;              % Echo time
SimPhantom.ppmaxis                  = SimParam.ppm_axis;
SimPhantom.Param.CSIdims            = SimPhantom.Param.dims(1:3);
SimPhantom.Param.NP                 = SimPhantom.Param.dims(end);
SimPhantom.Param.GyromagRatio       = SimParam.GyromagRatio;
SimPhantom.Param.B0                 = SimParam.B0;
SimPhantom.Param.ppmwindow          = SimParam.ppmwindow;
SimPhantom.Param.FirstOrdPhaseFunct = 0;
SimPhantom.HDOfreq                  = zeros(SimPhantom.Param.CSIdims) + 4.7;
SimPhantom.Param
   
%% Quantify with AMARES
SimPhantom.AMARESresults = ...
    AMARES4LowRank_v2(SimPhantom.img_t, SimPhantom.Param, SimPhantom.ppmaxis, SimPhantom.HDOfreq);

%% Generate parameter maps from AMARES fits
SimPhantom.AMARESmaps    = GenerateAMARESmaps(SimPhantom.AMARESresults);

end