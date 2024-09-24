function [fiddata_phacorr]=Phasecorrection(fiddata)
%% Phase correction using second point of the FID
dims=size(fiddata);
CSIdims=dims(2:end);
Ndata=prod(CSIdims);
NP=dims(1);
fidvec=reshape(fiddata,[NP Ndata]);

phase=atan2(imag(fidvec(2,:)),real(fidvec(2,:)));
fid_phcorr_vec=fidvec.*(exp(-1i*phase));

fiddata_phacorr=reshape(fid_phcorr_vec,dims);
end