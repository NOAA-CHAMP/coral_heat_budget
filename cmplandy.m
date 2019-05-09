function [stn,clim] = cmplandy(stn,PFX)
%function [stn,clim] = cmplandy(stn,PFX)
%
% Compare heat budget terms having prefix PFX (DEFAULT: 'ndbc_ncep_30a') with
% monthly CORE v.2 heat budget of Large and Yeager (2009).
%
% Last Saved Time-stamp: <Wed 2010-10-06 12:17:05  Lew.Gramer>

  if ( ~exist('PFX','var') || ~ischar(PFX) )
    PFX = 'ndbc_ncep_30a';
  end;
  if ( ~isempty(PFX) && PFX(end) ~= '_' )
    PFX(end+1) = '_';
  end;

  [stn,clim] = station_load_landy(stn);

  %DEBUG:  tic,

  [stn.monthly_sea_t,clim.sea_t,stn.monthly_sea_t_anom] = ...
      monthly_clim_ts(stn.ndbc_sea_t);

  srfld = 'ncep_srf';
  lrfld = 'ncep_lrf';
  lhfld = [PFX 'latent_heat_flux'];
  shfld = [PFX 'sensible_heat_flux'];
  [stn.monthly_sr,clim.sr,stn.monthly_sr_anom] = monthly_clim_ts(stn.(srfld));
  [stn.monthly_lr,clim.lr,stn.monthly_lr_anom] = monthly_clim_ts(stn.(lrfld));
  [stn.monthly_lh,clim.lh,stn.monthly_lh_anom] = monthly_clim_ts(stn.(lhfld));
  [stn.monthly_sh,clim.sh,stn.monthly_sh_anom] = monthly_clim_ts(stn.(shfld));

  %DEBUG:  toc,

return;
