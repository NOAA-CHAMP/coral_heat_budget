1;
% SCRIPT to analyze two-point gradients and three-point double-gradients in
% weekly AVHRR SST sea temperature fields from USF

ds=-3:3;

if ( ~isfield(stn,'avhrr_weekly_sst_field_rotated') )
  stn = station_rotate_field(stn,'isobath_orientation','avhrr_weekly_sst_field',[],true);
end;

clear x;

for ix=1:numel(ds)
  d = ds(ix);
  x.s(ix).lon=stn.avhrr_weekly_sst_field_rotated.lon(9,9+d);
  x.s(ix).lat=stn.avhrr_weekly_sst_field_rotated.lat(9,9+d);
  x.s(ix).date=stn.avhrr_weekly_sst_field_rotated.date;
  x.s(ix).data=stn.avhrr_weekly_sst_field_rotated.field(:,9,9+d);
end;

% Simple two-point gradient
for ix=1:numel(x.s)-1
  x.g(ix) = ts_op(x.s(ix+1),x.s(ix),'-');
  x.g(ix).data = x.g(ix).data ./ (stn.avhrr_weekly_sst_field_rotated.dx*1e3);
end;
% Simple three-point "Laplacian"
for ix=1:numel(x.g)-1
  x.l(ix) = ts_op(x.g(ix+1),x.g(ix),'-');
  x.l(ix).data = x.l(ix).data ./ (stn.avhrr_weekly_sst_field_rotated.dx*1e3);
end;


plot_ngdc_bathy_station(stn);
plot(stn.lon,stn.lat,'wp');
for ix=1:numel(ds)
  plot(x.s(ix).lon,x.s(ix).lat,'ws');
  text(x.s(ix).lon,x.s(ix).lat,['  \leftarrow  ',num2str(ds(ix)),'km'], 'Color','white');
end;
axis([stn.lon-0.1,stn.lon+0.1,stn.lat-0.1,stn.lat+0.1]); daspect([cosd(25),1,1]);



c='kbgrmcykbgrmcykbgrmcy';
s='.......sssssss^^^^^^^';

% accfun = @get_month;
accfun = @get_week;

% fmg; plot_ts(x.s); titlename([stn.station_name ' Temperatures']);
fmg;
for ix=1:numel(ds)
  [cum,tid] = grp_ts(x.s(ix).data,x.s(ix).date,accfun,@nanmean,0);
  plot(tid,cum,[c(ix) s(ix) '-']);
end;
legend(num2str(ds(:)'), 'Location','South');
titlename([stn.station_name ' Temperatures']);

%fmg; plot_ts(x.g); titlename([stn.station_name ' Gradients']);

fmg;
for ix=1:numel(x.g)
  [cum,tid] = grp_ts(x.g(ix).data,x.g(ix).date,accfun,@nanmean,0);
  plot(tid,cum,[c(ix) s(ix) '-']);
end;
legend(num2str(ds(1:end-1)'), 'Location','North');
titlename([stn.station_name ' Gradients']);

fmg;
for ix=1:numel(x.l)
  [cum,tid] = grp_ts(x.l(ix).data,x.l(ix).date,accfun,@nanmean,0);
  plot(tid,cum,[c(ix) s(ix) '-']);
end;
legend(num2str(ds(1:end-2)'), 'Location','South');
titlename([stn.station_name ' \partial_x_x']);


clear c cum d ds ix s tid x
