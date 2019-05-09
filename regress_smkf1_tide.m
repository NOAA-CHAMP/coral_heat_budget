1;

load('smkf1_tide_regression.mat');

[ix1,ix2]=intersect_dates(station.ndbc_tide_i_depth.date,station.tmd_tide.date);
x=station.ndbc_tide_m.data(ix1)-nanmean(station.ndbc_tide_m.data(ix1));
y=station.tmd_tide.data(ix2)-nanmean(station.tmd_tide.data(ix2));

[B,Stats] = robustfit(x,y);

fh = figure;
hold on;
plot(x,y,'b.');
plot(x,(B(1)+(B(2).*x)),'r-');
legs2 = sprintf('%.5g + %.5g*X, R^2~%0.2g, \n RMSE=%g, p=%0.5g, N=%d', ...
                  B(1), B(2), (Stats.coeffcorr(2,1)^2), ...
                  Stats.s, roundn(Stats.p(2),-5), length(x));
plot(x,x,'k-');
legend( 'NDBC vs. Mex', legs2, 'Perfect fit', 'Location','NorthWest');
