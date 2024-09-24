function ShowAmaresfit_v2(inputspec,DMIResults,row,col,slice,xaxis,Parameters,Referenceloc)

fitResults=DMIResults{row,col,slice};
time=Parameters.time+Parameters.TE;
dims=size(inputspec);
CSIdims=dims(2:end);
Freqshift=(Referenceloc(sub2ind(CSIdims, row, col, slice))-4.7);
for metabnum=1:numel(fitResults.amplitude)
    lineshapes(:,metabnum)=exp( -abs(time)*fitResults.linewidth(metabnum).' * pi).*exp(-pi^2*time.^2*(fitResults.sigma(metabnum).').^2);
    Resultfid(:,metabnum)=(fitResults.amplitude(metabnum).*lineshapes(:,metabnum).*exp(2*pi*1i*(fitResults.chemShift(metabnum))*(Parameters.Freq/(10^6))*time).'.*exp(0*1i*((fitResults.phase(metabnum)/180)*pi)));
end
figure('WindowState','maximized')
sig=inputspec(:,row,col,slice);%.*exp(-1i*((fitResults.phase(1)/180)*pi)).*Parameters.FirstOrdPhaseFunct;
fit=fftshift(fft(sum(Resultfid,2))).*Parameters.FirstOrdPhaseFunct  * 1 / sqrt(size(inputspec,1));

plot(xaxis-Freqshift,real(sig),'k','LineWidth',2)%Raw
hold on;
% plot(xaxis-Freqshift,imag(sig),'--k','LineWidth',2)%Imaginary Raw

plot(xaxis-Freqshift,real(fit),'b','LineWidth',2)%Fit
% plot(xaxis-Freqshift,imag(fit),'--b','LineWidth',2)%Imaginary Fit

plot(xaxis-Freqshift,real(sig-fit),'r','LineWidth',2)%Residual
% plot(xaxis-Freqshift,imag(sig-fit),'--r','LineWidth',2)%Imaginary Residual

plotshift1=max(real(inputspec(:,row,col,slice).*Parameters.FirstOrdPhaseFunct));
plotshift=plotshift1;

for mm=1:size(Resultfid,2)
    met=fftshift(fft(Resultfid(:,mm))).*Parameters.FirstOrdPhaseFunct;
    
    plot(xaxis,(real(met)-plotshift) * 1 / sqrt(size(inputspec,1)),'LineWidth',2)
    plotshift=plotshift+plotshift1;
end
hold off;
title(['AMARES Fit'],'FontSize',22)

xlim([-2 10])
set(gca,'XDir','reverse')
pbaspect([1 1 1])
xlabel('^2H frequency (ppm)','FontSize',22)
ax = gca;
ax.XAxis.FontSize=22;
ax.YAxis.FontSize=1;

if numel(fitResults.amplitude)==4
    legend('Raw','Fit','Residual',['Lipid: ',num2str(fitResults.amplitude(1),3)],['Glx: ',num2str(fitResults.amplitude(2),3)],['Glucose: ',num2str(fitResults.amplitude(3),3)],['HDO: ',num2str(fitResults.amplitude(4),3)],'FontSize',20,'Location','EastOutside')
elseif numel(fitResults.amplitude)==2
    legend('Raw','Fit','Residual',['Lipid: ',num2str(fitResults.amplitude(1),3)],['HDO: ',num2str(fitResults.amplitude(2),3)],'FontSize',20,'Location','EastOutside')
elseif numel(fitResults.amplitude)==3
    legend('Raw','Fit','Residual',['Lipid: ',num2str(fitResults.amplitude(1),3)],['Glucose: ',num2str(fitResults.amplitude(2),3)],['HDO: ',num2str(fitResults.amplitude(3),3)],'FontSize',20,'Location','EastOutside')
else
    legend('Raw','Fit','Residual','FontSize',20,'Location','EastOutside')
    
end
end