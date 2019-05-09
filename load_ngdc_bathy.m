function stn = load_ngdc_bathy(stn,doPlot)
%function stn = load_ngdc_bathy(stn,doPlot)
%
% Load NGDC 3-arcsecond (~92m) resolution Coastal Relief Model bathymetric
% data for the given STN struct or STNM reef monitoring station name. A file
% [STNM 'ngdc_bathy.dat'] must already exist in the datapath. If optional arg
% DOPLOT is TRUE, use CONTOURF (v.) to plot the loaded sea floor topography.
%
% Last Saved Time-stamp: <Fri 2010-10-08 20:47:48 Eastern Daylight Time gramer>
%
%error('This function now superceded by MATLAB/ECOFORECASTS/GET_NGDC_BATHY_STATION.M');

error('This function now superceded by MATLAB/ECOFORECASTS/GET_NGDC_BATHY_STATION.M');

  datapath = get_thesis_path('../data');

  if ( ~exist('doPlot','var') || isempty(doPlot) )
    doPlot = false;
  end;

  stnm = stn.station_name;

  fname = fullfile(datapath,[stnm '_ngdc_bathy.dat']);
  x = load(fname);
  lon = x(:,1);
  lat = x(:,2);
  dat = x(:,3);
  x = []; clear x;

  [LON,LAT] = meshgrid(unique(lon),unique(lat));
  stn.ngdc_92m_bathy.lon = LON';
  stn.ngdc_92m_bathy.lat = LAT';

  % nlon=length(unique(LON));
  % nlat=length(unique(LAT));
  % stn.ngdc_92m_bathy.field = reshape(dat',[nlon nlat]);
  stn.ngdc_92m_bathy.field = griddata(lon,lat,dat,LON',LAT');

  if ( doPlot )
    figure;
    contourf(stn.ngdc_92m_bathy.lon,stn.ngdc_92m_bathy.lat,stn.ngdc_92m_bathy.field);
    % set(gca,'ydir','rev');
    maxigraph;
    colorbar;
    titlename(['NDGC bathymetry surrounding ' stnm]);
  end;

return;
