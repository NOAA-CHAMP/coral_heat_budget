function stn = station_load_nocs(stn_or_stnm)
%function stn = station_load_nocs(stn_or_stnm)
%
% Extract time series of monthly sea-surface heat fluxes from the National
% Oceanographic Centre at Southhampton (NOCS) climatology for station STN.
% Fields added to STN include .NOCS_CLOUD_COVER, .NOCS_LATENT_HEAT_FLUX,
% .NOCS_SENSIBLE_HEAT_FLUX, .NOCS_AIR_T, .NOCS_SEA_T, .NOCS_SPECHUMID,
% .NOCS_BAROM, .NOCS_LRF (net longwave), .NOCS_DRF (net shortwave),
% .NOCS_WIND_SPEED, .NOCS_NET_HEAT_FLUX, and .NOCS_HEAT_FLUX_TERM. In
% addition, use STATION_HEAT_FLUX (v.) to calculate bulk heat fluxes from
% NOCS variables, including STN.nocs_relhumid, .nocs_bulk_latent_heat_flux,
% ...sensible_heat_flux,...wind_stress,...net_heat_flux,...heat_flux_term.
%
% All time series are *interpolated to hourly* using SPLINE. Original monthly
% time series are retained in separate fields, e.g., STN.MONTHLY_NOCS_AIR_T.
%
% Dataset description: http://www.noc.soton.ac.uk/noc_flux/noc2.php
% Data downloaded 2010 AUG 02 from: http://dss.ucar.edu/datasets/ds260.3/
% "Preliminary extensions" to years 2007-2009 are included in time series.
%
% Last Saved Time-stamp: <Wed 2012-03-28 13:00:55  Lew.Gramer>

  set_more off;

  datapath = get_thesis_path('../data');
  nocspath = fullfile(datapath,'nocs');

  ALLVARS = { ...
      'at',	'air_t' ; ...
      'cldc',	'cloud_cover' ; ...
      'lhf',	'latent_heat_flux' ; ...
      'lw',	'lrf' ; ...
      'qair',	'spechumid' ; ...
      'shf',	'sensible_heat_flux' ; ...
      'slp',	'barom' ; ...
      'sst',	'sea_t' ; ...
      'sw',	'srf' ; ...
      'wspd',	'wind_speed' ; ...
            };

  waswarned = false;

  if ( ischar(stn_or_stnm) )
    stnm = stn_or_stnm;
    stn.station_name = stnm;
  elseif ( isstruct(stn_or_stnm) )
    stn = stn_or_stnm;
    stnm = stn.station_name;
  end;

  result = [];

  matfname = fullfile(datapath,sprintf('%s_nocs.mat',lower(stnm)));

  if ( exist(matfname,'file') )

    disp(['Loading NOCS climatology from ' matfname]);
    load(matfname,'result');

  else

   if ( ~isfield(stn,'lon') || ~isfield(stn,'lat') )
    [stn.lon,stn.lat,stn.depth] = get_station_coords(stnm);
   end;

   result.station_name = stn.station_name;
   result.lon = stn.lon;
   result.lat = stn.lat;
   result.depth = stn.depth;

   disp(['Subsetting NOCS station climatology from raw files']);

   for varix = 1:size(ALLVARS,1)

    var = ALLVARS{varix,1};
    fld = ['nocs_' ALLVARS{varix,2}];
    mofld = ['monthly_' fld];

    result.(mofld).date = [];
    result.(mofld).data = [];
    result.(fld).date = [];
    result.(fld).data = [];

    for yr = 1987:2009

      if ( ~exist('latix','var') )
        latix = round(90+result.lat+0.5);
        lon = result.lon; lon(lon > 180) = lon - 360;
        lonix = round(180+(lon)+0.5);
      end;

      % ../data/nocs/nocs_v2_0_shf_2006.nc
      % NOTE: Files for 2007-2009 were renamed to remove '_prelim' from name 
      ncfname = fullfile(nocspath,sprintf('nocs_v2_0_%s_%04d.nc', var, yr));
      %DEBUG:      ncfname,
      nc = mDataset(ncfname);
      if ( isempty(nc) )
        warning('Unable to open NOCS file "%s"', ncfname);
        continue;
      end;
      dat = nc{var}(1:end,latix,lonix);
      if ( all(isnan(dat(:))) )
        if ( ~waswarned )
          warning('No data at grid-point! Doing 3x3 spatial averages.');
          waswarned = true;
        end;
        datmo = nc{var}(1:end,latix-1:latix+1,lonix-1:lonix+1);
        for mo = 1:12
          dat(mo) = nanmean(datmo(mo,:));
        end;
      end;
      close(nc); clear nc;

      result.(mofld).date(end+1:end+12,1) = datenum(yr,[1:12],1)';
      result.(mofld).data(end+1:end+12,1) = dat(:);

      dat = []; clear dat;
    end; %for yr

    % Turn nice compact monthly climatologies into big, smooth hourly ones
    [yr,mo,dy] = datevec(result.(mofld).date(end));
    dts = [result.(mofld).date(1):(1/24):datenum(yr,mo,31,23,59,59)]';
    %result.(fld).data = interp1(result.(mofld).date,result.(mofld).data,dts,'linear');
    result.(fld).data = interp1(result.(mofld).date,result.(mofld).data,dts,'spline',nan);
    result.(fld).date = dts;
    result.(fld).date(isnan(result.(fld).data)) = [];
    result.(fld).data(isnan(result.(fld).data)) = [];

   end; %for varix

   % Unit conversion [m/s] to kts
   % Convert wind speed to knots (because we'll convert it back in HFBULKTC!)
   result.monthly_nocs_wind_speed.data = mps2kts(result.monthly_nocs_wind_speed.data);
   result.nocs_wind_speed.data = mps2kts(result.nocs_wind_speed.data);
   % Unit conversion [g/kg] to [kg/kg]
   result.monthly_nocs_spechumid.data = result.monthly_nocs_spechumid.data ./ 1000;
   result.nocs_spechumid.data = result.nocs_spechumid.data ./ 1000;

% %%DEBUG ONLY
%    disp(['Interim save to ' matfname]);
%    save(matfname,'result');
%   end;
%   if (1)

   result = station_spechumid_to_relhumid(result,'monthly_nocs_air_t',...
                                          'monthly_nocs_spechumid',...
                                          'monthly_nocs_relhumid');
   result = station_spechumid_to_relhumid(result,'nocs_air_t',...
                                          'nocs_spechumid',...
                                          'nocs_relhumid');

   result.monthly_nocs_net_heat_flux.date = result.monthly_nocs_srf.date;
   result.monthly_nocs_net_heat_flux.data = ...
       result.monthly_nocs_srf.data + result.monthly_nocs_lrf.data + ...
       result.monthly_nocs_latent_heat_flux.data + result.monthly_nocs_sensible_heat_flux.data;

   result.nocs_net_heat_flux.date = result.nocs_srf.date;
   result.nocs_net_heat_flux.data = ...
       result.nocs_srf.data + result.nocs_lrf.data + ...
       result.nocs_latent_heat_flux.data + result.nocs_sensible_heat_flux.data;

   result = station_heat_flux(result,'nocs_wind_speed','nocs_air_t', ...
                              'nocs_relhumid','nocs_barom','nocs_sea_t', ...
                              'nocs_srf','nocs_lrf','nocs_bulk');

   hfld = 'tmd_tide_i_depth';
   if ( ~isfield(result,hfld) || isempty(result.(hfld).date) )
     hfld = 'tpxo_tide_i_depth';
     if ( ~isfield(result,hfld) || isempty(result.(hfld).date) )
       hfld = [];
     end;
     %DEBUG:     hfld,
   end;
   result = station_heat_flux_term(result,'nocs_net_heat_flux', ...
                                   'nocs_heat_flux_term', ...
                                   'nocs_sea_t',[],hfld);
   result = station_heat_flux_term(result,'nocs_bulk_net_heat_flux', ...
                                   'nocs_bulk_heat_flux_term', ...
                                   'nocs_sea_t',[],hfld);

   disp(['Saving NOCS climatology to ' matfname]);
   save(matfname,'result');

  end; % if exist matfname else


  % NOTE: REPLACE all NOCS fields in STN with those in RESULT
  flds = grepstruct(result,'nocs_');
  for ix = 1:length(flds)
    fld = flds{ix};
    stn.(fld) = result.(fld);
  end;

  set_more;

return;
