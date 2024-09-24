function [TransverseImage , CoronalImage]=GeneratesmoothImage(InputMatrix)

Dims=size(InputMatrix);
InputMatrix(isnan(InputMatrix))=0;

NP1=Dims(1);% AP dimension
NP2=Dims(2);% RL dimension
NP3=Dims(3);% FH dimension
factor=3; % Gaussian kernel size

if numel(Dims)==3
    disp('Input data is single time point.')
    NPtime=1;% Dynamics
elseif numel(Dims)==4
    disp('Input data is multiple time point.')
    NPtime=Dims(4);% Dynamics
    
end

%% Transverse images
p1 = 1:1:NP1*factor; p1 = p1 - NP1*factor/2;
p2 = 1:1:NP2*factor; p2 = p2 - NP2*factor/2;
[x,y] = ndgrid(p1,p2);
GaussianScale = 1/(2*factor*factor);
ConvKernelTra = exp(-GaussianScale*(x.^2 + y.^2));

for Transverseslice=1:NP3
    
    for time=1:NPtime
        TransverseImage(:,:,Transverseslice,time)=conv2(interp2(InputMatrix(:,:,Transverseslice,time),factor),ConvKernelTra,'same');
    end
end
%% Coronal images

p2 = 1:1:NP2*factor; p2 = p2 - NP2*factor/2;
p3 = 1:1:NP3*factor; p3 = p3 - NP3*factor/2;

[y,z] = ndgrid(p2,p3);
GaussianScale = 1/(2*factor*factor);
ConvKernelCor = exp(-GaussianScale*(y.^2 + z.^2));
for Coronalslice=1:NP1
    
    for time=1:NPtime
        CoronalImage(Coronalslice,:,:,time)=conv2(interp2(squeeze(InputMatrix(Coronalslice,:,end:-1:1,time)).',factor),ConvKernelCor,'same');% FH is flipped
    end
end
