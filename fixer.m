function stn = fixer(stn,PFX)
%error('RUN-ONCE fix-up code!')
error('RUN-ONCE fix-up code!')

  if ( ~exist('PFX','var') || isempty(PFX) )
    PFX = 'ndbc_bulk_rh';
  end;
  if ( PFX(end) ~= '_' )
    PFX(end+1) = '_';
  end;


  stnm = stn.station_name;
  Q0f = Q0factor(stn.ndbc_sea_t.data,[],2);


  if ( ~isempty(strmatch('ndbc_bulk_rh',PFX)) )
    [stn,Q0f] = ...
        station_heat_flux(stn,'ndbc_wind1_speed','ndbc_air_t','ndbc_relhumid',...
                          'ndbc_barom','ndbc_sea_t','ncep_srf','ndbc_bulk_rh_lrf',...
                          'ndbc_bulk_rh');
    stn.ndbc_bulk_rh_heat_flux_sum.date = stn.ndbc_bulk_rh_heat_flux_term.date;
    stn.ndbc_bulk_rh_heat_flux_sum.data = cumsum(stn.ndbc_bulk_rh_heat_flux_term.data);

  elseif ( ~isempty(strmatch('ndbc_bulk_26',PFX)) )
    [stn,Q0f] = ...
        station_heat_flux(stn,'ndbc_wind1_speed','ndbc_air_t','ndbc_relhumid',...
                          'ndbc_barom','ndbc_sea_t','ncep_srf','ndbc_bulk_rh_lrf',...
                          'ndbc_bulk_26','ncep_dsrf','ndbc_bulk_rh_dlrf');
    stn.ndbc_bulk_26_heat_flux_sum.date = stn.ndbc_bulk_26_heat_flux_term.date;
    stn.ndbc_bulk_26_heat_flux_sum.data = cumsum(stn.ndbc_bulk_26_heat_flux_term.data);

  elseif ( ~isempty(strmatch('ndbc_bulk_30',PFX)) )
    [stn,Q0f] = ...
        station_heat_flux(stn,'ndbc_wind1_speed','ndbc_air_t','ndbc_relhumid',...
                          'ndbc_barom','ndbc_sea_t','ncep_srf','ndbc_bulk_rh_lrf',...
                          'ndbc_bulk_30','ncep_dsrf','ndbc_bulk_rh_dlrf','ncep_precip');
    stn.ndbc_bulk_30_heat_flux_sum.date = stn.ndbc_bulk_30_heat_flux_term.date;
    stn.ndbc_bulk_30_heat_flux_sum.data = cumsum(stn.ndbc_bulk_30_heat_flux_term.data);
  end;


  hf = [PFX 'srf'];
  ft = [PFX 'srf_term'];
  if ( isfield(stn,hf) )
    stn.(ft).date = stn.(hf).date;
    stn.(ft).data = (stn.(hf).data ./ Q0f) .* (60*60);
  end;

  hf = [PFX 'lrf'];
  ft = [PFX 'lrf_term'];
  if ( isfield(stn,hf) )
    stn.(ft).date = stn.(hf).date;
    stn.(ft).data = (stn.(hf).data ./ Q0f) .* (60*60);
  end;

  hf = [PFX 'latent_heat_flux'];
  ft = [PFX 'latent_flux_term'];
  if ( isfield(stn,hf) )
    stn.(ft).date = stn.(hf).date;
    stn.(ft).data = (stn.(hf).data ./ Q0f) .* (60*60);
  end;

  hf = [PFX 'sensible_heat_flux'];
  ft = [PFX 'sensible_flux_term'];
  if ( isfield(stn,hf) )
    stn.(ft).date = stn.(hf).date;
    stn.(ft).data = (stn.(hf).data ./ Q0f) .* (60*60);
  end;

  hf = [PFX 'rain_heat_flux'];
  ft = [PFX 'rain_flux_term'];
  if ( isfield(stn,hf) )
    stn.(ft).date = stn.(hf).date;
    stn.(ft).data = (stn.(hf).data ./ Q0f) .* (60*60);
  end;



%   s = annocs(stnm);
%   for fld = {'nocs_latent_heat_flux','nocs_lrf','nocs_sensible_heat_flux', ...
%              'nocs_srf','nocs_net_heat_flux','nocs_heat_flux_term',...
%              'nocs_heat_flux_sum',...
%              'nocs_bulk_latent_heat_flux','nocs_bulk_sensible_heat_flux', ...
%              'nocs_bulk_net_heat_flux','nocs_bulk_heat_flux_term',...
%              'nocs_bulk_heat_flux_sum',}
%     stn.(fld{:}) = s.(fld{:});
%   end;
%   s = []; clear s;

return;
