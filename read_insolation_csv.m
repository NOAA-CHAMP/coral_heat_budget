function result = read_insolation_csv(fname, stanm)
%function result = read_insolation_csv(fname, stanm)
%
% Load CSV file 'fname' containing in/out shortwave (insolation) and longwave
% radiation, PAR, and other values estimated by NESDIS from satellite data,
% for gridpoints surrounding the reef monitoring site named 'STANM'. RESULT
% is empty, or is a struct with fields date, parW, par, su, sd, lu, ld, cld.
%
% Last Saved Time-stamp: <Mon 2011-10-24 16:03:10  lew.gramer>

  m = importdata(fname);
  dts = datenum(m.textdata(:,4));
  lats = str2num(strvcat(m.textdata(:,2)));
  lons = str2num(strvcat(m.textdata(:,3)));

  % dlon = min(diff(unique(lons)));
  % dlat = min(diff(unique(lats)));
  dlon = 0.125;
  dlat = 0.125;

  parW = m.data(:,2);
  par = m.data(:,3);
  sd = m.data(:,4);
  su = m.data(:,5);
  ld = m.data(:,6);
  lu = m.data(:,7);
  cld = m.data(:,8);

  m = []; clear m;

  %%%%%%% GOOD PLACE TO DO SOME GROSS QUALITY CONTROL!


  % Take a median of all nearby grid-points for each date-time
  [result.lon, result.lat] = get_station_coords(stanm);
  result.date = unique(dts);
  for dtix = 1:length(result.date)
    dt = result.date(dtix);
    ix = find( (dts == dt) & ...
               (abs(lons - result.lon) <= dlon) & ...
               (abs(lats - result.lat) <= dlat) );
    result.parW(dtix,1) = nanmedian(parW(ix));
    result.par(dtix,1) = nanmedian(par(ix));
    result.sd(dtix,1) = nanmedian(sd(ix));
    result.su(dtix,1) = nanmedian(su(ix));
    result.ld(dtix,1) = nanmedian(ld(ix));
    result.lu(dtix,1) = nanmedian(lu(ix));
    result.cld(dtix,1) = nanmedian(cld(ix));
  end;

return;
