% Save figures or output images for solved displacement and strain fields

%%
% Find img name
[~,imgname,imgext] = fileparts(file_name{1,ImgSeqNum});

%%
if DICpara.MethodToSaveFig == 1
    %% jpg
    figure(3-2); if DICpara.OrigDICImgTransparency == 0, colormap jet; caxis auto; end
    print([imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_DispU'],'-djpeg','-r300');
    
    figure(4-2); if DICpara.OrigDICImgTransparency == 0, colormap jet; caxis auto; end
    print([imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_DispV'],'-djpeg','-r300');
    
    figure(5-2); if DICpara.OrigDICImgTransparency == 0, colormap jet; caxis([-0.025,0.025]); end
    print([imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_strain_exx'],'-djpeg','-r300')
    
    figure(6-2); if DICpara.OrigDICImgTransparency == 0, colormap jet; caxis([-0.025,0.025]); end
    print([imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_strain_exy'],'-djpeg','-r300')
    
    figure(7-2); if DICpara.OrigDICImgTransparency == 0, colormap jet; caxis([-0.015,0.015]); end
    print([imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_strain_eyy'],'-djpeg','-r300')
    
    figure(8-2); if DICpara.OrigDICImgTransparency == 0, colormap jet;  caxis([0,0.025]); end
    print([imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_strain_principal_max'],'-djpeg','-r300')
    
    figure(9-2); if DICpara.OrigDICImgTransparency == 0, colormap jet;  caxis([-0.025,0]); end
    print([imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_strain_principal_min'],'-djpeg','-r300')
    
    figure(10-2); if DICpara.OrigDICImgTransparency == 0, colormap jet;  caxis([0,0.07]); end
    print([imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_strain_maxshear'],'-djpeg','-r300')
    
    figure(11-2); if DICpara.OrigDICImgTransparency == 0, colormap jet;  caxis([0,0.07]); end
    print([imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_strain_vonMises'],'-djpeg','-r300')
    
    
elseif DICpara.MethodToSaveFig == 2
    %% pdf
    filename = [imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_DispU'];
    figure(3-2); if DICpara.OrigDICImgTransparency == 0, colormap jet; caxis auto; end
    export_fig( gcf , '-pdf' , '-r300' , '-painters' , filename);
    
    filename = [imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_DispV'];
    figure(4-2); if DICpara.OrigDICImgTransparency == 0, colormap jet; caxis auto; end
    export_fig( gcf , '-pdf' , '-r300' , '-painters' , filename);
    
    filename = [imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_strain_exx'];
    figure(5-2); if DICpara.OrigDICImgTransparency == 0, colormap coolwarm(32); caxis([-0.025,0.025]); end
    export_fig( gcf , '-pdf' , '-r300' , '-painters' , filename);
    
    filename = [imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_strain_exy'];
    figure(6-2); if DICpara.OrigDICImgTransparency == 0, colormap coolwarm(32); caxis([-0.025,0.025]); end
    export_fig( gcf , '-pdf' , '-r300' , '-painters' , filename);
    
    filename = [imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_strain_eyy'];
    figure(7-2); if DICpara.OrigDICImgTransparency == 0, colormap coolwarm(32); caxis([-0.015,0.015]); end
    export_fig( gcf , '-pdf' , '-r300' , '-painters' , filename);
    
    filename = [imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_strain_principal_max'];
    figure(8-2); if DICpara.OrigDICImgTransparency == 0, colormap coolwarm(32); caxis([0,0.025]); end
    export_fig( gcf , '-pdf' , '-r300' , '-painters' , filename);
    
    filename = [imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_strain_principal_min'];
    figure(9-2); if DICpara.OrigDICImgTransparency == 0, colormap coolwarm(32); caxis([-0.025,0]); end
    export_fig( gcf , '-pdf' , '-r300' , '-painters' , filename);
    
    filename = [imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_strain_maxshear'];
    figure(10-2); if DICpara.OrigDICImgTransparency == 0, colormap coolwarm(32); caxis([0,0.07]); end
    export_fig( gcf , '-pdf' , '-r300' , '-painters' , filename);
    
    filename = [imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_strain_vonMises'];
    figure(11-2); if DICpara.OrigDICImgTransparency == 0, colormap coolwarm(32); caxis([0,0.07]); end
    export_fig( gcf , '-pdf' , '-r300' , '-painters' , filename);
    
    
else
    %% fig
    fprintf('Please modify codes manually in Section 8.');
    figure(3-2); colormap(coolwarm(128)); caxis auto; 
        savefig([imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_DispU.fig']);
    figure(4-2); colormap(coolwarm(128)); caxis auto; 
        savefig([imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_DispV.fig']);
    figure(5-2); colormap(coolwarm(128)); caxis([-0.05,0.1]); 
        savefig([imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_strain_exx.fig']);
    figure(6-2); colormap(coolwarm(128)); caxis([-0.05,0.05]); 
        savefig([imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_strain_exy.fig']);
    figure(7-2); colormap(coolwarm(128)); caxis([-0.1,0.05]); 
        savefig([imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_strain_eyy.fig']);
    figure(8-2); colormap(coolwarm(128)); caxis([0,0.025]); 
        savefig([imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_strain_principal_max.fig']);
    figure(9-2); colormap(coolwarm(128)); caxis([-0.025,0]); 
        savefig([imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_strain_principal_min.fig']);
    figure(10-2); colormap(coolwarm(128)); caxis([0,0.07]); 
        savefig([imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_strain_maxshear.fig']);
    figure(11-2); colormap(coolwarm(128)); caxis([0,0.07]); 
        savefig([imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_strain_vonMises.fig']);
    
    
end





