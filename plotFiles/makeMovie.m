files = dir('*_strain_eyy.jpg'); ImgGrayScaleMax = 255;
im = cell(length(files),1);
for i = 1:length(files)
    im{i} = files(i).name;
end

v = VideoWriter('video_strain_eyy.mp4');
v.FrameRate = 10;
open(v);
figure,
for tempk = [ 1:1:length(im) ]
    clf; 
    %imshow(imread(im{tempk},1),'DisplayRange',[0,ImgGrayScaleMax]);
    imshow(imread(im{tempk}) );  title(['Frame #',num2str(tempk+1)]); 
    % text(830,100,['Frame #',num2str(tempk+1)]);
    frame = getframe(gcf);
    writeVideo(v,frame);
    %waitbar(tempk/length(files));
    
end
close(v);

%%


alpha = 10:1:1000;
R = 66;
R0 = 25; 

F11 = alpha.^2./(alpha.^3+R^3-R0^3).^(2/3);
figure, plot(alpha,F11,'.'); 
