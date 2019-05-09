function frac = station_sand_cover(stn_or_stnm)
%function frac = station_sand_cover(stn_or_stnm)
%
% Estimate decimal fraction of benthic habitat which is sand or coral rubble (as
% opposed to live cover or hard-bottom), at the site named STN.station_name.
% Sources in literature include Moyer et al. 2003; Palandro et al. 2008; Walker
% et al. 2008; Bertelsen et al. 2009; Moses et al. 2009. These fractions are in
% turn used to estimate bottom reflectance, based on per-bottom-type reflectance
% data for coral reefs from Hochberg et al. (2003, 2004, pers. comm.)
%
% Last Saved Time-stamp: <Wed 2013-04-10 17:02:42 Eastern Daylight Time gramer>

  stn = get_station_from_station_name(stn_or_stnm);

  switch ( lower(stn.station_name) ),
   case 'lkwf1', frac = 0.40;
   case 'fwyf1', frac = 0.50;
   case 'ncorc', frac = 0.40;	% NCORE Site C
   case 'klgf1', frac = 0.40;	% NCORE Site Key Largo
   case 'mlrf1', frac = 0.50;
   case 'lonf1', frac = 0.30;
   case 'smkf1', frac = 0.50;
   case 'looe1', frac = 0.40;
   case 'sanf1', frac = 0.40;
   %case 'dryf1', frac = 0.30;
   case 'dryf1', frac = 0.50;
   otherwise,    error('Unrecognized station name "%s"!',stn.station_name);
  end;

return;
