1;

position = get(0, 'ScreenSize');
if ( position(4) >= 890 )
  % Lew's Desktop
  position = [1 61 1280 890];
else
  % Lew's Laptop
  position = [1 1 1280 726];
end;


imgfname = 'MODIS.2009003.160747.florida.chlor_a.png';
% imgfname = 'MODIS.2009003.161239.flbay.rgb.png';

rawimg = imread(imgfname);
filtimg = rawimg(650:900, 900:1200);
% filtimg = rawimg;
filtimg(filtimg > 20) = nan;
filtimg(filtimg < 10) = 0;
% filtimg = 255 - filtimg;

circen = houghcircle_range_of_radii(filtimg, 12, 30, 0.3, 36);

% [accum, circen, cirrad] = ...
%     CircularHough_Grd(filtimg, [12 36], 1, 18, 1.0);
% figure(1); imagesc(accum); axis image;
% title('Accumulation Array from Circular Hough Transform');
% set(gcf, 'Position', position);

figure(2); imagesc(filtimg); colormap('gray'); axis image;
hold on;
 plot(circen(:,1), circen(:,2), 'r+');
 for ix = 1:size(circen,1)
   %DrawCircle(circen(ix,1), circen(ix,2), cirrad(ix), 6, 'r');
   DrawCircle(circen(ix,1), circen(ix,2), circen(ix,3), 6, 'r');
 end;
hold off;
title('Raw Image with Circles Detected (center positions marked)');
set(gcf, 'Position', position);
