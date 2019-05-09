function stn = trim_erai_station(stn)
%function stn = trim_erai_station(stn)
%
% Remove the ERA Interim fields from station struct STN that are least
% likely to be useful for analysis: this reduces STN from 180 to 120Mb.
%
% Last Saved Time-stamp: <Wed 2011-01-19 20:54:22  lew.gramer>

  stn = rmfield(stn,grepstruct(stn,'erai_charnock'));
  stn = rmfield(stn,grepstruct(stn,'erai_barom'));
  stn = rmfield(stn,grepstruct(stn,'erai_conv_precip'));
  stn = rmfield(stn,grepstruct(stn,'erai_instant_wind_stress'));
  stn = rmfield(stn,grepstruct(stn,'erai_sunshine_duration'));
  stn = rmfield(stn,grepstruct(stn,'erai_clear_sky'));

return;
