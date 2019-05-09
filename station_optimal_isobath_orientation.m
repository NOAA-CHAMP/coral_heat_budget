function stn = station_optimal_isobath_orientation(stn_or_stnm)
%function stn = station_optimal_isobath_orientation(stn_or_stnm)
%
% Returns orientation of local isobath relative to True, estimated using
% a combination of local NGDC Coastal Relief Model (3-arcsec resolution)
% bathymetry and Google Earth imagery, with way points separated by ~1km, on
% 2011 Mar 18. Returns struct STN with field STN.isobath_orientation.
%
% Last Saved Time-stamp: <Fri 2016-07-15 14:43:09 Eastern Daylight Time gramer>

  stn = get_station_from_station_name(stn_or_stnm);

  orifld = 'isobath_orientation';

  if ( isfield(stn,orifld) )
    warning('STN.%s already present',orifld);

  else
    switch (lower(stn.station_name)),
     case 'lkwf1',	stn.(orifld)=  2;
     case 'fwyf1',	stn.(orifld)=  2;
     case 'mlrf1',	stn.(orifld)= 54;
     % case 'lonf1',	stn.(orifld)=200;  %? Estimate from multi-km-field NGDC bathymetry
     case 'lonf1',	stn.(orifld)=  0;  %? Estimate from *wide-field* NGDC bathymetry
     case 'tnrf1',	stn.(orifld)= 65;  %? GoogleEarth indicates convoluted isobaths
     case 'smkf1',	stn.(orifld)= 65;
     case 'looe1',	stn.(orifld)= 73; % Lee & Williams (1999) estimate: Changed GET_LOOE1_ADCP
     case 'sanf1',	stn.(orifld)= 82;  % Estimated in Google Earth 2012 Feb 21

     % % case 'dryf1',	stn.(orifld)=225;  % Estimated in Google Earth 2011 Jun 02
     % case 'dryf1',	stn.(orifld)= 70;  % Estimated from NGDC bathymetry same day
     case 'dryf1',	stn.(orifld)= 58;  % Estimated from multi-km-field* NGDC bathymetry

     % case 'plsf1',	stn.(orifld)=150;  % Estimated in Google Earth 2011 Jun 02
     case 'plsf1',	stn.(orifld)= 20;  % Estimated from NGDC bathymetry same day

     case 'ncorc',	stn.(orifld)= 42;  % NCORE C estimated from NGDC bathy. 2013 Apr 10
     case 'klgf1',	stn.(orifld)= 55;  % NCORE Key Largo from NGDC bathy. 2013 Apr 10

     case 'bnpin',	stn.(orifld)= 15;   % All from NGDC bathymetry 2011 Jun 23
     case 'bnpmi',	stn.(orifld)=330;
     case 'bnpon',	stn.(orifld)=345;
     case 'bnppa',	stn.(orifld)= 15;

     % Needs to be checked?
     case 'bnpnn',	stn.(orifld)=  0;

     % Need to be checked
     case 'cryf1',	stn.(orifld)= 30;
     case 'amsf1',	stn.(orifld)= 75;

     case '42003',	stn.(orifld)=  0;  % Meaningless for deep-water buoys

     % Lew.Gramer@noaa.gov 2016 Jul 15: Added based on A. Soloviev estimates
     case 'sfocb',     stn.(orifild)=  2;  % 'c-buoy',[26.0728,-80.0878, 20, 2]
     case 'sfoeb',     stn.(orifild)=  3;  % 'e-buoy',[26.0695,-80.0768, 50, 3]
     case 'sfone',     stn.(orifild)=  3;  % 'ne-buoy',[26.0695,-80.0768, 50, 3]
     case 'sfonw',     stn.(orifild)= 10;  % 'nw_w-btm',[26.0705,-80.0942, 11, 10]
     case 'sfose',     stn.(orifild)=  8;  % 'se-buoy',[26.1398,-80.0533, 145, 8]
     case 'sfosw',     stn.(orifild)=  1;  % 'sw-buoy',[26.0327,-80.0918, 20, 1]
     case 'sfopi',     stn.(orifild)=  4;  % 'pier-cc',[26.0577,-80.1090, 0, 4]

     otherwise,		error('Isobath orientation unknown for station %s',stn.station_name);
    end;
  end;

return;
