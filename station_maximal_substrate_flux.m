function stn = station_maximal_substrate_flux(stn,tfld,hfld,dTfld)
%function stn = station_maximal_substrate_flux(stn,tfld,hfld,dTfld)
%
% Calculate time series of heat flux into/out of the reef substrate that
% *would be required* to balance the heat budget in STN.(DTFLD), compared to
% observed variability in sea temperature STN.(TFLD) at depth STN.(HFLD).
%
% Last Saved Time-stamp: <Thu 2011-01-20 16:49:46 Eastern Standard Time gramer>

  if ( ~exist('tfld','var') || isempty(tfld) )
    tfld = 'ndbc_sea_t';
  end;
  if ( ~exist('hfld','var') || isempty(hfld) )
    hfld = 'tmd_tide_i_depth';
  end;
  if ( ~exist('dTfld','var') || isempty(dTfld) )
    dTfld = 'ndbc_erai_30a_ww3_fkeys_qe_dt';
  end;

  % Calculate hourly sea temperature change - but Mind the Gap!
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

  stn.maximal_substrate_flux_term.date = delt.date(deltix);
  stn.maximal_substrate_flux_term.data = stn.(dTfld).data(dTix) - delt.data(deltix);

  stn.maximal_substrate_flux.date = delt.date(deltix);
  stn.maximal_substrate_flux.data = stn.maximal_substrate_flux_term.data .* Q0fac;

return;
