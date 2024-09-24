function ResultMaps=GenerateAMARESmaps(AMARESoutput)

AMARESmatrix=cell2mat(AMARESoutput);
Amp_map=zeros(numel(AMARESmatrix),size(AMARESmatrix(1).amplitude,2));
LW_map=zeros(numel(AMARESmatrix),size(AMARESmatrix(1).amplitude,2));
ChemShift_map=zeros(numel(AMARESmatrix),size(AMARESmatrix(1).amplitude,2));
tic
for m=1:numel(AMARESmatrix)
    Amp_map(m,:)=AMARESmatrix(m).amplitude;
    LW_map(m,:)=AMARESmatrix(m).linewidth;
    ChemShift_map(m,:)=AMARESmatrix(m).chemShift;
end
toc
ResultMaps.Amplitude=reshape(Amp_map,[size(AMARESmatrix) size(AMARESmatrix(1).amplitude,2)]);
ResultMaps.Linewidth=reshape(LW_map,[size(AMARESmatrix) size(AMARESmatrix(1).amplitude,2)]);
ResultMaps.ChemicalShift=reshape(ChemShift_map,[size(AMARESmatrix) size(AMARESmatrix(1).amplitude,2)]);
end
