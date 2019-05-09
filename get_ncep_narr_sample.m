1;

fname = '../data/ncep_narr_2008.mat';

if ( exist(fname,'file') )
  load(fname, 'dat2008');
else
  dat2008 = query_ncep_narr_subset([], datenum(2008,01,01), datenum(2008,12,31));
  disp(['Saving dat2008 to ' fname]);
  save(fname, 'dat2008');
end;

% Subset our world list of stations to those inside our BBOX
stns = get_all_station_metadata;
XV = [dat2008.bbox(3) dat2008.bbox(4) dat2008.bbox(4) dat2008.bbox(3)];
YV = [dat2008.bbox(2) dat2008.bbox(2) dat2008.bbox(1) dat2008.bbox(1)];
goodix = find( inside(stns.lons, stns.lats, XV, YV) );
clear XV YV;
stns.codes  = stns.codes(goodix);
stns.lons   = stns.lons(goodix);
stns.lats   = stns.lats(goodix);
stns.depths = stns.depths(goodix);


for ix = 1:length(stns.lons)

  stnm = lower(stns.codes{ix});
  disp(stnm);
  [lonix,latix] = gridnbhd_km(dat2008.lon,dat2008.lat,stns.lons(ix),stns.lats(ix),0);
  dat2008.(stnm).lonix = lonix;
  dat2008.(stnm).latix = latix;

  dat2008.(stnm).relhumid.date = dat2008.date;
  dat2008.(stnm).ncep_dsrf.date = dat2008.date;
  dat2008.(stnm).ncep_usrf.date = dat2008.date;
  dat2008.(stnm).ncep_dlrf.date = dat2008.date;
  dat2008.(stnm).ncep_ulrf.date = dat2008.date;
  for dix = 1:length(dat2008.date)
    dat2008.(stnm).relhumid.data = squeeze(dat2008.Relative_humidity(:,lonix,latix));
    dat2008.(stnm).ncep_dsrf.data = squeeze(dat2008.Downward_shortwave_radiation_flux(:,lonix,latix));
    dat2008.(stnm).ncep_usrf.data = squeeze(dat2008.Upward_short_wave_radiation_flux_surface(:,lonix,latix));
    dat2008.(stnm).ncep_dlrf.data = squeeze(dat2008.Downward_longwave_radiation_flux(:,lonix,latix));
    dat2008.(stnm).ncep_ulrf.data = squeeze(dat2008.Upward_long_wave_radiation_flux_surface(:,lonix,latix));
  end;

end;

disp(['Saving again to ' fname]);
save(fname, 'dat2008');
