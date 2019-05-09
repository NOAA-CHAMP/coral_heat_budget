1;

more off;

datapath = get_thesis_path('../data');

for yr = [ 2009 2007 2006 2005 2004 2003 2002 2001 2000 1999 1998 1997 1996 1995 1994 1993 1992 1991 1990 1989 1988 1987 ]

 switch ( yr )
  case 2009,
   mos = [1 2];
  case 2010,
   mos = [];
  otherwise,
   mos = 1:12;
 end;

 for mo = mos(:)'

  dat = []; clear dat;

  fname = fullfile(datapath, sprintf('ncep_cfsr_%04d_%02d.mat', yr, mo));

  disp(['Creating ' fname '...']);

  %DEBUG:    dat = query_ncep_cfsr_subset([], datenum(yr,mo,1), datenum(yr,mo,2));
  %%%% ??? DEBUG
  dat = query_ncep_cfsr_subset(yr, mo);
  disp(['Saving dat to ' fname]);
  save(fname, 'dat');

  % Subset our world list of stations to those inside our BBOX
  if ( ~exist('stns','var') || isempty(stns) )
    stns = get_all_station_metadata;
    XV = [dat.bbox(3) dat.bbox(4) dat.bbox(4) dat.bbox(3)];
    YV = [dat.bbox(2) dat.bbox(2) dat.bbox(1) dat.bbox(1)];
    goodix = find( inside(stns.lons, stns.lats, XV, YV) );
    clear XV YV;
    stns.codes  = stns.codes(goodix);
    stns.lons   = stns.lons(goodix);
    stns.lats   = stns.lats(goodix);
    stns.depths = stns.depths(goodix);
  end;


  for ix = 1:length(stns.lons)

    stnm = lower(stns.codes{ix});
    disp(stnm);
    [lonix,latix] = gridnbhd_km(dat.lon,dat.lat,stns.lons(ix),stns.lats(ix),0);
    dat.(stnm).lonix = lonix;
    dat.(stnm).latix = latix;

%     dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'cfsr_tau_u', 'wndstrs_u');
%     dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'cfsr_tau_v', 'wndstrs_v');
%     dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'cfsr_wind_u', 'wnd10m_u');
%     dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'cfsr_wind_v', 'wnd10m_v');
    dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'cfsr_usrf', 'uswsfc');
    dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'cfsr_ulrf', 'ulwsfc');
%     dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'cfsr_sea_t', 'tmpsfc');
%     dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'cfsr_air_t', 'tmp2m');
    dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'cfsr_sensible_heat_flux', 'shtfl');
    dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'cfsr_spechumid', 'q2m');
%     dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'cfsr_barom', 'pressfc');
    dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'cfsr_precip', 'prate');
%     dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'cfsr_ocean_surfcur_u', 'ocnv5');
%     dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'cfsr_ocean_surfcur_v', 'ocnu5');
%     dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'cfsr_ocean_ml_w', 'ocnvv55');
%     dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'cfsr_ocean_sst', 'ocnsst');
%     dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'cfsr_ocean_sss', 'ocnsal5');
%     dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'cfsr_ocean_mld', 'ocnmld');
    dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'cfsr_latent_heat_flux', 'lhtfl');
    dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'cfsr_dsrf', 'dswsfc');
    dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'cfsr_dlrf', 'dlwsfc');

  end;

  %%%% ??? DEBUG
  disp(['Saving again to ' fname]);
  save(fname, 'dat');

 end; %for mo

end; %for yr
