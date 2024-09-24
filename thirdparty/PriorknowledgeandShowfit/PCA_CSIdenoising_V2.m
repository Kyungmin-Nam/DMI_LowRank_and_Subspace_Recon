function denoisedCSI_spec=PCA_CSIdenoising_V2(spec_data,patch_size,Parameters)
%% PCA based denoising of CSI data
% Eliminate eigenvalues fitting marchenko pasteur distrubition

% Calculate eigenvalues for each patch. Eigen values coming from noise should represent a
% Marchenko-Pastur distribution. Calculate threshold and set eigen values
% below this threshold to zero. Same approach Martijn Froeling used.
disp('Starting PCA based denoising.')
disp(strcat('Patch size:',num2str(patch_size)))
tic
NP=size(spec_data,1);
extended_data=repmat(spec_data,[1 3 3 3]);
denoisedCSI_spec=spec_data;
n=patch_size^3;
% Denoise by elimination eigenvalues equal/below noise eigenvalues
for mm=1:prod(Parameters.CSIdims)
    [Grid_row,Grid_col,Grid_slice]=ind2sub([Parameters.CSIdims(1) Parameters.CSIdims(2) Parameters.CSIdims(3)],mm);
    signalpatch=zeros(patch_size^3,NP*2);
    
    for m=1:n
        [row,col,slice]=ind2sub([patch_size patch_size patch_size],m);
        signalpatch(m,1:NP)=real(extended_data(:,Grid_row+row+(Parameters.CSIdims(1)-ceil(patch_size/2)),Grid_col+col+(Parameters.CSIdims(2)-ceil(patch_size/2)),Grid_slice+slice+(Parameters.CSIdims(3)-ceil(patch_size/2))));
        signalpatch(m,NP+1:end)=imag(extended_data(:,Grid_row+row+(Parameters.CSIdims(1)-ceil(patch_size/2)),Grid_col+col+(Parameters.CSIdims(2)-ceil(patch_size/2)),Grid_slice+slice+(Parameters.CSIdims(3)-ceil(patch_size/2))));
    end
    denoisedpatch=DenoiseMatrix(signalpatch).';
    
    denoisedsignal=denoisedpatch(1:NP,ceil(n/2))+1i*denoisedpatch(NP+1:end,ceil(n/2));
    denoisedCSI_spec(:,Grid_row,Grid_col,Grid_slice)=denoisedsignal;
    
    clear denoisedsignal signalpatch U S V;
    
end

toc
disp('Finished PCA based denoising.')

end
