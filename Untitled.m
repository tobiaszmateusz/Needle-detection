% %  EXAMPLE #1:
% %  rawimg = imread('TestHT_Pentagon.png');
%  rawimg = dicomread('IM_0011.dcm','frames',100)
%  rawimg = rawimg(150:650,1:550)
% %  rawimg = imrotate(rawimg,14,'crop');
%  fltr4img = [1 2 3 2 1; 2 3 4 3 2; 3 4 6 4 3; 2 3 4 3 2; 1 2 3 2 1];
%  fltr4img = fltr4img / sum(fltr4img(:));
%  imgfltrd = filter2( fltr4img , rawimg );
%  tic;
%  [accum, axis_rho, axis_theta, lineprm, lineseg] = ...
%      Hough_Grd(imgfltrd, 6, 0.02);
%  toc;
%  figure(1); imagesc(axis_theta*(180/pi), axis_rho, accum); axis xy;
%  xlabel('Theta (degree)'); ylabel('Pho (pixels)');
%  title('Accumulation Array from Hough Transform');
%  figure(2); imagesc(rawimg); colormap('gray'); axis image;
%  DrawLines_2Ends(lineseg);
%  title('Raw Image with Line Segments Detected');
%  %%
%  %  EXAMPLE #2:
%  rawimg = dicomread('IM_0011.dcm','frames',100)
%  rawimg = rawimg(150:650,1:550)
%  rawimg = imrotate(rawimg,14,'crop');
%  fltr4img = [1 1 1 1 1; 1 2 2 2 1; 1 2 4 2 1; 1 2 2 2 1; 1 1 1 1 1];
%  fltr4img = fltr4img / sum(fltr4img(:));
%  imgfltrd = filter2( fltr4img , rawimg );
%  tic;
%  [accum, axis_rho, axis_theta, lineprm, lineseg] = ...
%      Hough_Grd(imgfltrd, 8, 0.05);
%  toc;
%  figure(1); imagesc(axis_theta*(180/pi), axis_rho, accum); axis xy;
%  xlabel('Theta (degree)'); ylabel('Pho (pixels)');
%  title('Accumulation Array from Hough Transform');
%  figure(2); imagesc(imgfltrd); colormap('gray'); axis image;
%  DrawLines_Polar(size(imgfltrd), lineprm);
%  title('Raw Image (Blurred) with Lines Detected');
%  figure(3); imagesc(rawimg); colormap('gray'); axis image;
%  DrawLines_2Ends(lineseg);
%  title('Raw Image with Line Segments Detected');
%  
%  %%
%  rawimg = dicomread('IM_0011.dcm','frames',100)
%  rawimg = rawimg(150:650,1:550,:)
% %  imshow(rawimg);
%  gray=rgb2gray(rawimg)
%  bw=imbinarize(gray)
%  bw=edge(bw)
%  imshow(bw)
%  [h,t,r]=hough(bw,'Theta',-90:1:89,'RhoResolution',1);
%  P=houghpeaks(h,6,'threshold',0.50,'NHoodSize',[1 7]);
%  lines=houghlines(bw,t,r,P,'FillGap',10,'MinLength',20);
%  imshow(bw)
%  hold on
%  for ii = 1:length(lines)
%      xy = [lines(ii).point1; lines(ii).point2];
%      plot(xy(:,1),xy(:,2),"LineWidth",2,"Color","green");
%  end
%  
%  



 %%
img = imread("C:\Users\Mateusz\Desktop\AOM\Dane do projektu\DICOM\test2\Series_001/IM_0009_Frame150.jpg");
 %%imshow(img)
togray = rgb2gray(img);
[J,rect] = imcrop(togray);
rect=[122.5100  280.5100  365.9800  296.9800];
croppedImage = imcrop(togray, rect);
imshow(croppedImage);
shock = shock_filter(croppedImage,11,35,0.28);


figure(2)
imshow(img)
figure(2)
imshow(shock)
%%
% [G,GABOUT]=gaborfilter(shock,0.3,0.15,0,0)
% 
% R=real(GABOUT);
% I=imag(GABOUT);
% M=abs(GABOUT);
% P=angle(GABOUT);
% 
% % Show the filter's outputs
% figure;
% colormap(redgreen);
% subplot(2,2,1);
% k=127.5/max(max(abs(R)));
% image(uint8(k*R+127.5));
% subplot(2,2,2);
% k=127.5/max(max(abs(I)));
% image(uint8(k*I+127.5));
% 
% 
% % Show the magnitudes
% figure;
% colormap(grayscale);
% k=255/max(max(M));
% image(uint8(k*M));
% 
% % Show the phases
% figure;
% colormap(redgreen);
% k=127.5/pi;
% image(uint8(k*P+127.5));
% 

%%
lambda  = 5;
theta   = 45;
psi     = [0 pi/2];
gamma   = 12;
bw      = 5;
N       = 16;

img_out = zeros(size(shock,1), size(shock,2), N);
for n=1:N
    gb = gabor_fn(bw,gamma,psi(1),lambda,theta)...
        + 1i * gabor_fn(bw,gamma,psi(2),lambda,theta);
    % gb is the n-th gabor filter
    img_out(:,:,n) = imfilter(shock, gb, 'symmetric');
    % filter output to the n-th channel
    theta = theta + 2*pi/N;
    % next orientation
end
figure(1);
imshow(shock);
title('input image');
figure(2);
img_out_disp = sum(abs(img_out).^2, 3).^0.5;
% default superposition method, L2-norm
img_out_disp = img_out_disp./max(img_out_disp(:));
% normalize
imshow(img_out_disp);
title('gabor output, L-2 super-imposed, normalized');

%%

%%
wynikprogowania = proguj(img_out_disp,0.5);
CC = bwconncomp(wynikprogowania,4)
numOfPixels = cellfun(@numel,CC.PixelIdxList)
[unused,indexOfMax] = max(numOfPixels)
biggest = zeros(size(wynikprogowania))
biggest(CC.PixelIdxList{indexOfMax}) = 1


CMap = [0,0,0; 1,0,0];
RGB  = ind2rgb(biggest + 1, CMap)
imshow(RGB)
%%
rgbImage = croppedImage;
[rows, columns, numberOfColorChannels] = size(rgbImage)
% Display the test image.
subplot(2, 2, 1);
imshow(rgbImage, []);
axis('on', 'image');
drawnow;
hp = impixelinfo(); % Set up status line to see values when you mouse over the image.
% Set up figure properties:
% Enlarge figure to full screen.
hFig1 = gcf;
hFig1.Units = 'Normalized';
hFig1.WindowState = 'maximized';
% Get rid of tool bar and pulldown menus that are along top of figure.
% set(gcf, 'Toolbar', 'none', 'Menu', 'none');
% Give a name to the title bar.
hFig1.Name = 'Demo by Image Analyst';

mask = RGB;
% Display the initial mask image.
subplot(2, 2, 2);
imshow(mask, []);
hp = impixelinfo(); % Set up status line to see values when you mouse over the image.
axis('on', 'image');
drawnow;
% Mask the image using bsxfun() function to multiply the mask by each channel individually.  Works for gray scale as well as RGB Color images.
maskedRgbImage = bsxfun(@times, rgbImage, cast(mask, 'like', rgbImage));
% Display the final masked image.
subplot(2, 2, 3);
imshow(maskedRgbImage, []);
axis('on', 'image');
drawnow;
hp = impixelinfo(); % Set up status line to see values when you mouse over the image.
% Display the final masked image of the background by inverting the mask.
% Mask the image using bsxfun() function to multiply the mask by each channel individually.  Works for gray scale as well as RGB Color images.
backgroundImage = bsxfun(@times, rgbImage, cast(~mask, 'like', rgbImage));
subplot(2, 2, 4);
imshow(backgroundImage, []);
axis('on', 'image');
drawnow;
hp = impixelinfo(); % Set up status line to see values when you mouse over the image.
%%
% imageSize = size(shock);
% numRows = imageSize(1);
% numCols = imageSize(2);
% 
% wavelengthMin = 4/sqrt(2);
% wavelengthMax = hypot(numRows,numCols);
% n = floor(log2(wavelengthMax/wavelengthMin));
% wavelength = 2.^(0:(n-2)) * wavelengthMin;
% 
% deltaTheta = 45;
% orientation = 0:deltaTheta:(180-deltaTheta);
% 
% g = gabor(wavelength,orientation);
% 
% %%
% gabormag = imgaborfilt(shock,g);
% 
% %%
% for i = 1:length(g)
%     sigma = 0.5*g(i).Wavelength;
%     K = 3;
%     gabormag(:,:,i) = imgaussfilt(gabormag(:,:,i),K*sigma); 
% end
% 
% %%
% X = 1:numCols;
% Y = 1:numRows;
% [X,Y] = meshgrid(X,Y);
% featureSet = cat(3,gabormag,X);
% featureSet = cat(3,featureSet,Y);
% 
% %%
% numPoints = numRows*numCols;
% X = reshape(featureSet,numRows*numCols,[]);
% 
% %%
% X = bsxfun(@minus, X, mean(X));
% X = bsxfun(@rdivide,X,std(X));
% 
% %%
% coeff = pca(X);
% feature2DImage = reshape(X*coeff(:,1),numRows,numCols);
% figure
% imshow(feature2DImage,[])
% 
% %%
% L = kmeans(X,2,'Replicates',5);
% 
% %%
% L = reshape(L,[numRows numCols]);
% figure
% imshow(label2rgb(L))