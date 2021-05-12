%%
   info = dicominfo('C:\Users\Mateusz\Desktop\AOM\Dane do projektu\DICOM\IM_0009.dcm');
   saveshock = ones(297, 366, info.NumberOfFrames);
   savecropped = zeros(297, 366, 3, info.NumberOfFrames, 'uint8');
   shock2 = ones(297, 366, info.NumberOfFrames);
   img_out_disp2 = ones(297, 366, info.NumberOfFrames);
   image_data = dicomread('C:\Users\Mateusz\Desktop\AOM\Dane do projektu\DICOM\IM_0009.dcm');
   rect=[122.5100  280.5100  365.9800  296.9800];
%    implay(image_data)
test = image_data(:,:,150);
for i = 1:info.NumberOfFrames
    rawimg = dicomread('IM_0009.dcm','frames',i);
    croppedImage2 = imcrop(rawimg, rect);
    togray2 = rgb2gray(croppedImage2);
    shock2 = shock_filter(togray2,11,35,0.28); 
    saveshock(:,:,i) = shock2;
    savecropped(:,:,:,i) = croppedImage2;
end

%%
lambda  = 5;
theta   = 45;
psi     = [0 pi/2];
gamma   = 12;
bw      = 5;
N       = 16;
for i = 1:info.NumberOfFrames
img_out = zeros(size(saveshock(:,:,i),1), size(saveshock(:,:,i),2), N);
for n=1:N
    gb = gabor_fn(bw,gamma,psi(1),lambda,theta)...
        + 1i * gabor_fn(bw,gamma,psi(2),lambda,theta);
    % gb is the n-th gabor filter
    img_out(:,:,n) = imfilter(saveshock(:,:,i), gb, 'symmetric');
    % filter output to the n-th channel
    theta = theta + 2*pi/N;
    % next orientation
    img_out_disp = sum(abs(img_out).^2, 3).^0.5;
    img_out_disp = img_out_disp./max(img_out_disp(:));
    img_out_disp2(:,:,i) = img_out_disp;
end
end

%%
biggest2 = ones(297, 366, info.NumberOfFrames);
rgb2 = ones(297, 366, 3, 255);
%%
% wynikprogowania = proguj(img_out_disp,0.5);
% CC = bwconncomp(wynikprogowania,4)
% numOfPixels = cellfun(@numel,CC.PixelIdxList)
% [unused,indexOfMax] = max(numOfPixels)
% biggest = zeros(size(wynikprogowania))
% biggest(CC.PixelIdxList{indexOfMax}) = 1


for i = 1:info.NumberOfFrames
    wynikprogowania = proguj(img_out_disp2(:,:,i),0.5);
    CC = bwconncomp(wynikprogowania,4);
    numOfPixels = cellfun(@numel,CC.PixelIdxList);
    [unused,indexOfMax] = max(numOfPixels);
    biggest = zeros(size(wynikprogowania));
    biggest(CC.PixelIdxList{indexOfMax}) = 1;
    biggest2(:,:,i) = biggest;
end

%%
for i = 1:info.NumberOfFrames
    CMap = [0,0,0; 1,0,0];
    RGB  = ind2rgb(biggest2(:,:,i) + 1, CMap);
    
    rgb2(:,:,:,i) = RGB;
end


%%
backgroundImage = zeros(297, 366, 3, info.NumberOfFrames, 'uint8');
backgroundImage = zeros(297, 366, 3, info.NumberOfFrames, 'uint8');
%%
for i = 1:info.NumberOfFrames
mask = rgb2(:,:,:,i);

maskedRgbImage(:,:,:,i) = bsxfun(@times, savecropped(:,:,:,i), cast(mask, 'like', savecropped(:,:,:,i)));

backgroundImage(:,:,:,i) = bsxfun(@times, savecropped(:,:,:,i), cast(~mask, 'like', savecropped(:,:,:,i)));
end

%%

implay(backgroundImage);

%%
% wspi = zeros(608,1);
% wspj = zeros(608,1);
% temp = 1;
% for i = 1:297
%     for j = 1:366
%         if(biggest2(i,j,150)==1)
%             wspi(temp) = i;
%             wspj(temp) = j;
%             temp = temp + 1;
%         end
%     end
% end
% wsp1a = wspi(1,1);
% wsp1b = wspi(608,1);
% 
% wsp2a = wspj(1,1);
% wsp2b = wspj(608,1);
% 
% imshow(biggest2(:,:,150))
% hold on
% plot(wsp1b,wsp1a,"ro")
% plot(wsp2b,wsp2a,"ro")

%%
wspi = zeros(608,1);
wspja = zeros(6,1);
wspjb = zeros(6,1);
temp = 1;
for i = 1:297
    for j = 1:366
        if(biggest2(i,j,150)==1)
            wspi(temp) = i;
            temp = temp + 1;
        end
    end
end

%%
temp = 1;
for j = 1:366
    if(biggest2(wspi(2,1),j,150)==1)
        wspja(temp) = j;
        temp = temp + 1;
    end
end
temp = 1;
for j = 1:366
	if(biggest2(wspi(607,1),j,150)==1)
        wspjb(temp) = j;
        temp = temp + 1;
	end
end
%%
wsp1a = wspi(1,1);
wsp1b = wspi(608,1);

wsp2a = wspja(1,1)-2;
wsp2b = round(mean(wspjb));

figure()
imshow(biggest2(:,:,150))
hold on
plot(wsp1b,wsp1a,"ro")
plot(wsp2b,wsp2a,"ro")
%%
a = (wsp1a-wsp2a)/(wsp1b-wsp2b);
b = wsp1a - a*wsp1b;

x=0:1:366;
y=a*x+b;

figure()
imshow(img_out_disp2(:,:,150))
hold on
plot(x,y,"r");
figure()
imshow(backgroundImage(:,:,:,150))
hold on
plot(x,y,"r");

%%
% % HoG detector
% [featureVector,hogVisualization] = extractHOGFeatures(img_out_disp2(:,:,150),'CellSize',[16 16]);
% figure;
% imshow(img_out_disp2(:,:,150)); 
% hold on;
% plot(hogVisualization);

%%
cellSize = 16;
hog1 = vl_hog(single(img_out_disp2(:,:,1)),cellSize, 'verbose', 'variant', 'dalaltriggs');
imhog1 = vl_hog('render', hog1, 'verbose', 'variant', 'dalaltriggs');
clf ; imagesc(imhog1) ; colormap gray ;

%%
cellSize = 16;
hog2 = vl_hog(single(img_out_disp2(:,:,150)),cellSize, 'verbose', 'variant', 'dalaltriggs');
imhog2 = vl_hog('render', hog2, 'verbose', 'variant', 'dalaltriggs');
clf ; imagesc(imhog2) ; colormap gray ;

%%
hogfinal = hog1 - hog2;