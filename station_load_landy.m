function [stn,clim] = station_load_landy(stn_or_stnm)
%function [stn,clim] = station_load_landy(stn_or_stnm)
%
% Load flux terms (short-/longwave radiative, sensible/latent turbulent,
% tau_x/y wind stress, runoff/evaporation/precipitation) from monthly CORE
% v.2 heat budget of Large and Yeager (2009). Return station struct STN with
% '.landy_*' fields added, and CLIM struct with just '.landy_*' fields.
%
% Last Saved Time-stamp: <Sat 2010-12-11 17:46:46 Eastern Standard Time gramer>

  datapath = get_thesis_path('../data');
  corepath = fullfile(datapath,'CORE');

  if ( ischar(stn_or_stnm) )
    stn.station_name = stn_or_stnm;
  elseif ( isfield(stn_or_stnm,'station_name') )
    stn = stn_or_stnm;
  else
    error('First arg must be station name or struct with field .station_name!');
  end;

  if ( ~isfield(stn,'lat') )
    [stn.lon,stn.lat,stn.depth] = get_station_coords(stn.station_name);
  end;

  %DEBUG:  tic,

  clim = [];

  matfname = fullfile(datapath,sprintf('%s_landy.mat',lower(stn.station_name)));

  if ( exist(matfname,'file') )
    disp(sprintf('Loading station Large&Yeager climatology from "%s"',matfname));
    load(matfname,'clim');

  else
    disp('Subsetting station Large&Yeager climatology from raw netCDF...');

    yrs = [1987:2006];

    for yr = yrs(:)'

      ncfname = fullfile(corepath,sprintf('CORE.v2.%04d.nc',yr));
      nc = mDataset(ncfname);
      %DEBUG:  nj_info(nc),
      lats = nc{'lat'}(:,:);
      lons = nc{'lon'}(:,:);
      % Idiotic negative West longitudes
      lons(lons>180) = lons(lons>180) - 360;
      [ig,lonix]=min(abs(stn.lon-lons));
      [ig,latix]=min(abs(stn.lat-lats));

      % Latent, longwave (D and U), sensible, shortwave surface fluxes, resp.

      clim = append_landy(clim,yr,'landy_lh',(nc{'Q_lat'}(:,latix,lonix)));
      clim = append_landy(clim,yr,'landy_lrd',(nc{'Q_lwdn'}(:,latix,lonix)));
      % NCEP (and other sources) provide "POSITIVE" upward longwave flux
      clim = append_landy(clim,yr,'landy_lru',(- nc{'Q_lwup'}(:,latix,lonix)));
      clim = append_landy(clim,yr,'landy_sh',(nc{'Q_sen'}(:,latix,lonix)));
      clim = append_landy(clim,yr,'landy_sr',(nc{'Q_swnet'}(:,latix,lonix)));
      clim = append_landy(clim,yr,'landy_ev',(nc{'F_evap'}(:,latix,lonix)));
      clim = append_landy(clim,yr,'landy_pr',(nc{'F_prec'}(:,latix,lonix)));
      % Runoff is in [kg/s/m^2] - like evaporation and precipitation
      clim = append_landy(clim,yr,'landy_ro',(nc{'F_roff'}(:,latix,lonix)));
      clim = append_landy(clim,yr,'landy_taux',(nc{'taux'}(:,latix,lonix)));
      clim = append_landy(clim,yr,'landy_tauy',(nc{'tauy'}(:,latix,lonix)));
      close(nc); clear nc;

    end;

    clim.landy_lr.date = clim.landy_lrd.date;
    clim.landy_lr.data = clim.landy_lrd.data - clim.landy_lru.data;

    clim.landy_net_heat_flux.date = clim.landy_lh.date;
    clim.landy_net_heat_flux.data = clim.landy_lh.data + clim.landy_sh.data + clim.landy_sr.data + clim.landy_lr.data;

    disp(sprintf('Saving station Large&Yeager climatology to "%s"',matfname));
    save(matfname,'clim');

  end;

  flds = fieldnames(clim);
  for fldix = 1:length(flds)
    fld = flds{fldix};
    stn.(fld) = clim.(fld);
  end;

  %DEBUG:  toc,

return;

function clim = append_landy(clim,yr,fldnm,dat)
  if ( ~isfield(clim,fldnm) )
    clim.(fldnm).date = [];
    clim.(fldnm).data = [];
  end;
  dts = datenum(yr,[1:12],1);
  clim.(fldnm).date(end+1:end+12,1) = dts(:);
  clim.(fldnm).data(end+1:end+12,1) = dat(:);
return;
