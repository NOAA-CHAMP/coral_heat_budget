function stn = trim_looe1_adcp(stn)
%function stn = trim_looe1_adcp(stn)
%
% Remove the ADCP fields from Looe-Key station struct STN that are least
% likely to be useful for analysis: this reduces STN from 400 to 230Mb!
%
% Last Saved Time-stamp: <Wed 2011-01-19 17:30:14 Eastern Standard Time gramer>

  stn = rmfield(stn,grepstruct(stn,'adcp_baroclinic'));
  stn = rmfield(stn,grepstruct(stn,'adcp_err'));
  stn = rmfield(stn,grepstruct(stn,'adcp_pct_good'));
  stn = rmfield(stn,grepstruct(stn,'eacnt'));
  stn = rmfield(stn,grepstruct(stn,'adcp_w'));

return;
