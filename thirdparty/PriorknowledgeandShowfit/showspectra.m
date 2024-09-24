function showspectra(SPECTRA,ppm_axis)

specmax=max(real(SPECTRA),[],'all');
zlimitwindow=find(ppm_axis>0 & ppm_axis<10);%

specmin=min(real(SPECTRA(zlimitwindow,:,:)),[],'all');
nkx=size(SPECTRA,2);
nky=size(SPECTRA,3);

figure('color','w','WindowState','maximized')
for idx_iy = 1:nky
    for idx_ix = 1:nkx
        %                         im_spectra = double(real(SPECTRA(:,idx_ix, idx_iy).*exp(-1i * (2* pi * (ppm_axis-4.7).'*45.7 * 0.00184))));% With first order phase correction
        %         im_spectra = double(real(SPECTRAStraight(:,idx_ix, idx_iy)));
        %         im_spectra2 = double(real(SPECTRAbended(:,idx_ix, idx_iy)));
        im_spectra = double(real(SPECTRA(:,idx_ix, idx_iy)));
        
        hAxe = axes(...
            'Parent'    , gcf                        , ...
            'Box'       , 'on'                       , ...
            'LineWidth' , 4                          , ...
            'Position'  , [(idx_iy-1)*1/nky, 1-(idx_ix*1/nkx), 1/nky, 1/nkx], ...
            'XDir'      , 'reverse'                  , ...
            'XTick'     , []                         , ...%'XTick'     , [0 2 4 6 8 10]
            'YTick'     , []                         , ...
            'XTickLabel', {''; ''}                   , ...
            'YTickLabel', {''; ''}                   , ...
            'TickDir'   , 'out'                      , ...
            'Ylim'      , [specmin*1.2 specmax*1.2]     , ...
            'XLim'      , [0 10]                     ,...
            'XColor'    , [1 0 0]                    ,...
            'YColor'    , [1 0 0]                    ,...
            'Color'     , [0 0 0]);% for black background [0 0 0] ---- for white background [1 1 1] 
        %         if idx_iy*idx_ix==4
        %             hLine = line( ...
        %                 'Parent', hAxe            , ...
        %                 'XData' , ppm_axis-0.3        , ...
        %                 'YData' , im_spectra      , ...
        %                 'Linewidth', 3         , ...
        %                 'Color',[0 0 0]);
        %             %         hold on
        %             %         hLine = line( ...
        %             %             'Parent', hAxe            , ...
        %             %             'XData' , ppm_axis        , ...
        %             %             'YData' , im_spectra2      , ...
        %             %             'Linewidth', 2         , ...
        %             %             'Color',[1 0 0]);
        %             %         hold off
        %             % %         Add line on 4.7 ppm
        %         else
        hLine = line( ...
            'Parent', hAxe            , ...
            'XData' , ppm_axis        , ...
            'YData' , im_spectra      , ...
            'Linewidth', 3.2         , ...
            'Color',[0 1 0]);% for black lines [0 0 0] ---- for green lines [0 1 0]
        %         end
        waterline=0;
        if waterline==1
                hold on
                WaterLine = line( ...
                    'Parent', hAxe            , ...
                    'XData' , [4.7 4.7]        , ...
                    'YData' , [specmin*1.5 specmax*1.2]      , ...
                    'Linewidth', 4         , ...
                    'LineStyle', '--' , ...
                    'Color',[0 0 1]);
                hold off
        end
    end
end
end