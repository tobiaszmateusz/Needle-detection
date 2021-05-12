[X, map] = dicomread('IM_0011.dcm');
% mov = immovie(X,map)
% implay(mov)

X2 = dicomread('IM_0011.dcm','frames',100)
X3 = X2(:,1:600)

X3 = imrotate(X3,14,'crop');
X4 = edge(X3,'canny');
figure(1)
imshow(X3)
figure(2)
imshow(X4)
%%
[H,theta,rho] = hough(X4);


figure(3)
imshow(imadjust(rescale(H)),[],...
       'XData',theta,...
       'YData',rho,...
       'InitialMagnification','fit');
xlabel('\theta (degrees)')
ylabel('\rho')
axis on
axis normal 
hold on
colormap(gca,hot)

P = houghpeaks(H,11,'threshold',ceil(0.3*max(H(:))));


x = theta(P(:,2));
y = rho(P(:,1));
plot(x,y,'s','color','black');


lines = houghlines(X4,theta,rho,P,'FillGap',13,'MinLength',9);


figure, imshow(X3), hold on
max_len = 0;
for k = 1:length(lines)
   xy = [lines(k).point1; lines(k).point2];
   plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

   % Plot beginnings and ends of lines
   plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
   plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');

   % Determine the endpoints of the longest line segment
   len = norm(lines(k).point1 - lines(k).point2);
   if ( len > max_len)
      max_len = len;
      xy_long = xy;
   end
end
% highlight the longest line segment
plot(xy_long(:,1),xy_long(:,2),'LineWidth',2,'Color','red');

%%
[X, map] = dicomread('IM_0011.dcm');
I = dicomread('IM_0011.dcm','frames',100)
I = I(:,1:600)

BW = edge(I,'canny');
[H,T,R] = hough(BW,'RhoResolution',0.5,'Theta',-90:0.5:89);

subplot(2,1,1);
imshow(I);
title('gantrycrane.png');
subplot(2,1,2);
imshow(imadjust(rescale(H)),'XData',T,'YData',R,...
      'InitialMagnification','fit');
title('Hough transform of gantrycrane.png');
xlabel('\theta'), ylabel('\rho');
axis on, axis normal, hold on;
colormap(gca,hot);


P = houghpeaks(H,11,'threshold',ceil(0.3*max(H(:))));


x = T(P(:,2));
y = R(P(:,1));
plot(x,y,'s','color','black');


lines = houghlines(X4,T,R,P,'FillGap',5,'MinLength',7);


figure, imshow(X3), hold on
max_len = 0;
for k = 1:length(lines)
   xy = [lines(k).point1; lines(k).point2];
   plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');

   % Plot beginnings and ends of lines
   plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
   plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');

   % Determine the endpoints of the longest line segment
   len = norm(lines(k).point1 - lines(k).point2);
   if ( len > max_len)
      max_len = len;
      xy_long = xy;
   end
end
% highlight the longest line segment
plot(xy_long(:,1),xy_long(:,2),'LineWidth',2,'Color','red');