function stn = station_substrate_flux(stn,tfld,hfld,dTfld)
%function stn = station_substrate_flux(stn,tfld,hfld,dTfld)
%
% Calculate time series of heat flux into/out of the reef substrate.
%
% Last Saved Time-stamp: <Mon 2011-01-24 17:34:29  lew.gramer>
%error('Function replaced by changes to STATION_ABSORBED_INSOLATION!');

error('Function replaced by changes to STATION_ABSORBED_INSOLATION!');

  if ( ~exist('tfld','var') || isempty(tfld) )
    tfld = 'ndbc_sea_t';
  end;
  if ( ~exist('hfld','var') || isempty(hfld) )
    hfld = 'tmd_tide_i_depth';
  end;
  if ( ~exist('dTfld','var') || isempty(dTfld) )
    dTfld = 'ndbc_erai_30a_ww3_fkeys_qe_dt';
  end;

  % Ab: Sea-bottom mean reflectance
  % (Based on mix of sea-bottom habitat types)
  try
    frac = station_sand_cover(stn);
  catch
    % Hard-bottom/sea grass mean >=30%, living reef down to <10%
    frac = 0.6;
  end;
  Ab = (0.40 * frac) + (0.07 * (1-frac));

  Qswi = ???


  delt.date = stn.(tfld).date(1:end-1);
  delt.data = diff(stn.(tfld).data);
  gapix = find(diff(delt.date) > (1.5/24));
  delt.date(gapix) = [];
  delt.data(gapix) = [];

  [deltix,tix,hix,dTix] = ...
      intersect_all_dates([],delt.date,stn.(tfld).date,stn.(hfld).date,stn.(dTfld).date);

  % Cp and rho from [Fofonoff and Millard, 1983]
  s = repmat(36, size(delt.data(deltix)));
  t = stn.(tfld).data(tix);
  h = stn.(hfld).data(hix);
  rho = sw_dens0( s, t );
  Cp = sw_cp( s, t, h );
  Q0fac = (rho.*Cp.*h)./3600;

  stn.substrate_flux_term.date = delt.date(deltix);
  stn.substrate_flux_term.data = stn.(dTfld).data(dTix) - delt.data(deltix);

  stn.substrate_flux.date = delt.date(deltix);
  stn.substrate_flux.data = stn.substrate_flux_term.data .* Q0fac;

return;
