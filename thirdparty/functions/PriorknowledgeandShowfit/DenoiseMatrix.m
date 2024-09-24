function OutMatrix=DenoiseMatrix(InpMatrix)
%% Debug
% InpMatrix=noisepatch;
% InpMatrix=signalpatch;
%
n=size(InpMatrix,1);
m=size(InpMatrix,2);
r=min(m,n);
% mean_data = mean(InpMatrix); InpMatrix=InpMatrix-mean_data;
[U,S,V] = svd(InpMatrix,'econ');
eigv = flip((diag(S)).^2);
lam_r = eigv(1)/ n ;
clam = 0;
sigma2 = NaN;
tolerance=1;
for p=1:r % Grid search
    lam = eigv(p) / n;
    clam = clam+lam;
    gam = (m-p)/n;
    sigsq1 = (clam/(p));
    sigsq2 = (lam-lam_r)/(4*sqrt(gam));
    if(sigsq2 < sigsq1)
        sigma2 = sqrt(sigsq1+sigsq2)/2;% In case if you want to generate a noise map
        cutoff_p = p-tolerance;
    end
end

cutoff_p = r-cutoff_p;
eigv = flip(eigv);

if(cutoff_p > 1)
    Snew = zeros(size(S));
    Sd = diag(sqrt(eigv(1:cutoff_p)));
    Snew(1:cutoff_p,1:cutoff_p) = Sd;
    rebuilt_data = reshape(U*Snew*V',size(InpMatrix));%+mean_data;
    
else
    rebuilt_data=InpMatrix;
end
% noise=sigma2

% showS_Eig=1; Turn it on for Debug
% if showS_Eig==1
%     figure
%     subplot(2,3,1)
%     histogram(S)
%     title('Histogram of input S')
%     subplot(2,3,4)
%     plot(diag(S))
%     hold on
%     if(cutoff_p > 1)
%         plot(diag(Snew))
%     end
%     hold off
%     title('plot of S')
%     
%     subplot(1,3,2)
%     histogram(eigv,500)
%     title('Histogram of input eigenvalues')
%     
%     subplot(1,3,3)
%     histogram(eigv(1:cutoff_p),500)
%     title('Histogram of output eigenvalues')
% end
OutMatrix=rebuilt_data;
