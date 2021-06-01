% Save figures or output images

%%
% Find img name
[~,imgname,imgext] = fileparts(file_name{1,ImgSeqNum});

%%
if DICpara.MethodToSaveFig == 1
    %% jpg
    figure(3); if DICpara.OrigDICImgTransparency == 0, colormap jet; caxis auto; end
    caxis auto; print([imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_DispU'],'-djpeg','-r300');
    
    figure(4); if DICpara.OrigDICImgTransparency == 0, colormap jet; caxis auto; end
    caxis auto; print([imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_DispV'],'-djpeg','-r300');
    
    figure(5); if DICpara.OrigDICImgTransparency == 0, colormap jet; caxis([-0.025,0.025]); end
    print([imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_exx'],'-djpeg','-r300')
    
    figure(6); if DICpara.OrigDICImgTransparency == 0, colormap jet; caxis([-0.025,0.025]); end
    print([imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_exy'],'-djpeg','-r300')
    
    figure(7); if DICpara.OrigDICImgTransparency == 0, colormap jet; caxis([-0.015,0.015]); end
    print([imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_eyy'],'-djpeg','-r300')
    
    figure(8); if DICpara.OrigDICImgTransparency == 0, colormap jet;  caxis([0,0.025]); end
    print([imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_principalstrain_max'],'-djpeg','-r300')
    
    figure(9); if DICpara.OrigDICImgTransparency == 0, colormap jet;  caxis([-0.025,0]); end
    print([imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_principalstrain_min'],'-djpeg','-r300')
    
    figure(10); if DICpara.OrigDICImgTransparency == 0, colormap jet;  caxis([0,0.07]); end
    print([imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_maxshear'],'-djpeg','-r300')
    
    figure(11); if DICpara.OrigDICImgTransparency == 0, colormap jet;  caxis([0,0.07]); end
    print([imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_vonMises_strain'],'-djpeg','-r300')
    
    
elseif DICpara.MethodToSaveFig == 2
    %% pdf
    filename = [imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_DispU'];
    figure(3); if DICpara.OrigDICImgTransparency == 0, colormap jet; caxis auto; end
    export_fig( gcf , '-pdf' , '-r300' , '-painters' , filename);
    
    filename = [imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_DispV'];
    figure(4); if DICpara.OrigDICImgTransparency == 0, colormap jet; caxis auto; end
    export_fig( gcf , '-pdf' , '-r300' , '-painters' , filename);
    
    filename = [imgname,'_WS',num2str(winsize),'_ST',num2str(DICpara.winstepsize),'_exx'];
    figure(5); if DICpara.OrigDICImgTransparency == 0, colormap coolwarm(32); caxis([-0.025,0.025]); end
    export_fig( gcf , '-pdf' , '-r300' , '-painters' , filename);
    
    filename = [imgname,'_WS',num2str(winsize),'_ST',num2str(DICpara.winstepsize),'_exy'];
    figure(6); if DICpara.OrigDICImgTransparency == 0, colormap coolwarm(32); caxis([-0.025,0.025]); end
    export_fig( gcf , '-pdf' , '-r300' , '-painters' , filename);
    
    filename = [imgname,'_WS',num2str(winsize),'_ST',num2str(DICpara.winstepsize),'_eyy'];
    figure(7); if DICpara.OrigDICImgTransparency == 0, colormap coolwarm(32); caxis([-0.015,0.015]); end
    export_fig( gcf , '-pdf' , '-r300' , '-painters' , filename);
    
    filename = [imgname,'_WS',num2str(winsize),'_ST',num2str(DICpara.winstepsize),'_principalstrain_max'];
    figure(8); if DICpara.OrigDICImgTransparency == 0, colormap coolwarm(32); caxis([0,0.025]); end
    export_fig( gcf , '-pdf' , '-r300' , '-painters' , filename);
    
    filename = [imgname,'_WS',num2str(winsize),'_ST',num2str(DICpara.winstepsize),'_principalstrain_min'];
    figure(9); if DICpara.OrigDICImgTransparency == 0, colormap coolwarm(32); caxis([-0.025,0]); end
    export_fig( gcf , '-pdf' , '-r300' , '-painters' , filename);
    
    filename = [imgname,'_WS',num2str(winsize),'_ST',num2str(DICpara.winstepsize),'_maxshear'];
    figure(10); if DICpara.OrigDICImgTransparency == 0, colormap coolwarm(32); caxis([0,0.07]); end
    export_fig( gcf , '-pdf' , '-r300' , '-painters' , filename);
    
    filename = [imgname,'_WS',num2str(winsize),'_ST',num2str(DICpara.winstepsize),'_vonMises_strain'];
    figure(11); if DICpara.OrigDICImgTransparency == 0, colormap coolwarm(32); caxis([0,0.07]); end
    export_fig( gcf , '-pdf' , '-r300' , '-painters' , filename);
    
    
else
    %% fig
    fprintf('Please modify codes manually in Section 8.');
    figure(3); colormap(coolwarm(128)); caxis auto; savefig([imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_DispU.fig']);
    figure(4); colormap(coolwarm(128)); caxis auto; savefig([imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_DispV.fig']);
    figure(5); colormap(coolwarm(128)); caxis([-0.05,0.1]); savefig([imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_exx.fig']);
    figure(6); colormap(coolwarm(128)); caxis([-0.05,0.05]); savefig([imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_exy.fig']);
    figure(7); colormap(coolwarm(128)); caxis([-0.1,0.05]); savefig([imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_eyy.fig']);
    figure(8); colormap(coolwarm(128)); caxis([0,0.025]); savefig([imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_principalstrain_max.fig']);
    figure(9); colormap(coolwarm(128)); caxis([-0.025,0]); savefig([imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_principalstrain_min.fig']);
    figure(10); colormap(coolwarm(128)); caxis([0,0.07]); savefig([imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_maxshear.fig']);
    figure(11); colormap(coolwarm(128)); caxis([0,0.07]); savefig([imgname,'_WS',num2str(DICpara.winsize),'_ST',num2str(DICpara.winstepsize),'_vonMises_strain.fig']);
    
    
end





