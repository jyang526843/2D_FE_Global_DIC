files = dir('*_DispV.jpg'); ImgGrayScaleMax = 255;
im = cell(length(files),1);
for i = 1:length(files)
    im{i} = files(i).name;
end

v = VideoWriter('frac_heter_dispv.mp4');
open(v);
figure,
for tempk = 1 : length(files)
    clf
    %imshow(imread(im{tempk},1),'DisplayRange',[0,ImgGrayScaleMax]);
    imshow(imread(im{tempk}) ); 
   % hold on;
   % viscircles([CircleFitPar(tempk,2),CircleFitPar(tempk,1)],Rnew(tempk),'Color','b');
    %axis([1,512,1,128]);
    frame = getframe(gcf);
    writeVideo(v,frame);
    %waitbar(tempk/length(files));
    
end
close(v);