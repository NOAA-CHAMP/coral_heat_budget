function stn = station_calc_quasi_eulerian(stn,stokesu,stokesv,meanu,meanv,eulerpfx)
%function stn = station_calc_quasi_eulerian(stn,stokesu,stokesv,meanu,meanv,eulerpfx)
%
% Combine Stokes drift and mean current data into quasi-Eulerian current
% components U, V, speed, and direction. EULERPFX is name prefix for quasi-
% Eulerian currents (DEFAULT 'quasi_eulerian'). New fields added to STN are
% [EULERPFX '_u'], [EULERPFX '_v'], [EULERPFX '_speed'], [EULERPFX '_dir'].
%
% Last Saved Time-stamp: <Wed 2011-05-04 08:35:00  lew.gramer>

  if ( ~exist('eulerpfx','var') || isempty(eulerpfx) )
    eulerpfx = 'quasi_eulerian';
  end;

  euleru = [eulerpfx '_u'];
  eulerv = [eulerpfx '_v'];
  eulers = [eulerpfx '_speed'];
  eulerd = [eulerpfx '_dir'];

  [hix,wix] = intersect_dates(stn.(meanu).date,stn.(stokesu).date);

  stn.(euleru).date = stn.(meanu).date(hix);
  stn.(euleru).data = stn.(meanu).data(hix) + stn.(stokesu).data(wix);
  stn.(eulerv).date = stn.(meanv).date(hix);
  stn.(eulerv).data = stn.(meanv).data(hix) + stn.(stokesv).data(wix);
  stn.(eulers).date = stn.(euleru).date;
  stn.(eulers).data = uv_to_spd(stn.(euleru).data,stn.(eulerv).data);
  stn.(eulerd).date = stn.(euleru).date;
  stn.(eulerd).data = uv_to_dir_curr(stn.(euleru).data,stn.(eulerv).data);

return;
