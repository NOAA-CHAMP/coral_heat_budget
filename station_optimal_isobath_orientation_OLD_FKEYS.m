function stn = station_optimal_isobath_orientation(stn)
%function stn = station_optimal_isobath_orientation(stn)
%
% Returns orientation of local isobath relative to True, which optimizes the
% heat advection term as calculated from FKEYS HYCOM model outputs 2004-2008.
% Struct STN is returned with field STN.isobath_orientation populated. NOTE
% these orientations were also confirmed using Google Earth imagery with way
% points separated by approximately 1km, on 2011 Mar 18.
%
% Last Saved Time-stamp: <Tue 2012-02-21 13:57:27  lew.gramer>

  orifld = 'isobath_orientation';
  if ( ~isfield(stn,orifld) )
    switch (lower(stn.station_name)),
     case 'fwyf1',	stn.(orifld)=1.80;
     % case 'mlrf1',	stn.(orifld)=52.25;
     case 'mlrf1',	stn.(orifld)=54.20;  %INTERP_FKEYS
     case 'lonf1',	stn.(orifld)=67.50;
     case 'tnrf1',	stn.(orifld)=65;  %? GoogleEarth indicates convoluted isobaths
     case 'smkf1',	stn.(orifld)=65.00;
     case 'looe1',	stn.(orifld)=77;  % Straight isobaths in GoogleEarth: Changed GET_LOOE1_ADCP
     otherwise,		error('Isobath orientation unknown for station %s',stn.station_name);
    end;
  end;

return;
