function station = calc_quasi_eulerian(station,stokespfx,meanpfx,eulerpfx)
%function station = calc_quasi_eulerian(station,stokespfx,meanpfx,eulerpfx)
%
% Combine Stokes drift and mean current data to construct quasi-Eulerian
% current components U, V, speed, and direction. Last three args are string
% prefixes for each of Stokes, mean, and quasi-Eulerian current fields.
%
% Last Saved Time-stamp: <Sun 2010-07-25 17:44:31 Eastern Daylight Time gramer>


  if ( ~exist('stokespfx','var') || isempty(stokespfx) )
    stokespfx = 'stokes_drift';
  end;
  if ( ~exist('meanpfx','var') || isempty(meanpfx) )
    meanpfx = 'global_hycom';
  end;
  if ( ~exist('eulerpfx','var') || isempty(eulerpfx) )
    eulerpfx = 'quasi_eulerian';
  end;

  meanu = [meanpfx '_u'];
  meanv = [meanpfx '_v'];
  stokesu = [stokespfx '_u'];
  stokesv = [stokespfx '_v'];

  euleru = [eulerpfx '_u'];
  eulerv = [eulerpfx '_v'];
  eulers = [eulerpfx '_speed'];
  eulerd = [eulerpfx '_dir'];

  [hix,wix] = intersect_dates(station.(meanu).date,station.(stokesu).date);

  station.(euleru).date = station.(meanu).date(hix);
  station.(euleru).data = station.(meanu).data(hix) + station.(stokesu).data(wix);
  station.(eulerv).date = station.(meanv).date(hix);
  station.(eulerv).data = station.(meanv).data(hix) + station.(stokesv).data(wix);
  station.(eulers).date = station.(euleru).date;
  station.(eulers).data = uv_to_spd(station.(euleru).data,station.(eulerv).data);
  station.(eulerd).date = station.(euleru).date;
  station.(eulerd).data = uv_to_dir_curr(station.(euleru).data,station.(eulerv).data);

return;
