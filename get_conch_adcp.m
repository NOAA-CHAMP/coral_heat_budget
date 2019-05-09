function stn = get_conch_adcp(stn)
%function stn = get_conch_adcp(stn)
%
% Return all data from SFP Conch Reef ADCP deployment in station struct STN.
%
% Last Saved Time-stamp: <Sun 2013-02-24 14:30:38 Eastern Standard Time gramer>

  if ( ~exist('stn','var') || isempty(stn) )
    stn = get_station_from_station_name('aqua1');
  end;

  cnch = load(fullfile(get_thesis_path('../data'),'SFP','conch_hourly_adcp.mat'));

  dts = datenum([cnch.adcp.year],[cnch.adcp.month],[cnch.adcp.day],[cnch.adcp.hour],[cnch.adcp.minute],[cnch.adcp.second])';

  stn.adcp_i_depth.date = dts;
  stn.adcp_i_depth.data = [cnch.adcp.water_depth]';


  vav = [cnch.adcp.velocity_avg];
  vmn = [vav.vertical_mean];
  sfc = [vav.sfc3bins];
  btm = [vav.btm3bins];
  vav=[]; clear vav;

  vbn = [cnch.adcp.velocity_bin];

  wvs = [cnch.adcp.waves];
  cnch=[]; clear cnch;


  stn.adcp_speed.date = dts;
  stn.adcp_speed.data = [vmn.mag]'./100;
  stn.adcp_dir.date = dts;
  stn.adcp_dir.data = [vmn.dir]';
  stn.adcp_u.date = dts;
  stn.adcp_u.data = [vmn.u]'./100;
  stn.adcp_v.date = dts;
  stn.adcp_v.data = [vmn.v]'./100;
  vmn=[]; clear vmn;

  stn.adcp_sfc_speed.date = dts;
  stn.adcp_sfc_speed.data = [sfc.mag]'./100;
  stn.adcp_sfc_dir.date = dts;
  stn.adcp_sfc_dir.data = [sfc.dir]';
  stn.adcp_sfc_u.date = dts;
  stn.adcp_sfc_u.data = [sfc.u]'./100;
  stn.adcp_sfc_v.date = dts;
  stn.adcp_sfc_v.data = [sfc.v]'./100;
  sfc=[]; clear sfc;

  stn.adcp_btm_speed.date = dts;
  stn.adcp_btm_speed.data = [btm.mag]'./100;
  stn.adcp_btm_dir.date = dts;
  stn.adcp_btm_dir.data = [btm.dir]';
  stn.adcp_btm_u.date = dts;
  stn.adcp_btm_u.data = [btm.u]'./100;
  stn.adcp_btm_v.date = dts;
  stn.adcp_btm_v.data = [btm.v]'./100;
  btm=[]; clear btm;

  n = numel(vbn)/50;
  if ( n ~= floor(n) )
    error('Shape error for currents profiles');
  end;
  emix=find(arrayfun(@(x)(isempty(x.mag)|isempty(x.dir)|isempty(x.u)|isempty(x.v)),vbn));
  for emixix=emix(:)';
    vbn(emixix).mag = nan;
    vbn(emixix).dir = nan;
    vbn(emixix).u = nan;
    vbn(emixix).v = nan;
  end;

  stn.adcp_speed.prof = reshape([vbn.mag]',[50,n])';
  stn.adcp_dir.prof = reshape([vbn.dir]',[50,n])';
  stn.adcp_u.prof = reshape([vbn.u]',[50,n])';
  stn.adcp_v.prof = reshape([vbn.v]',[50,n])';
  vbn=[]; clear vbn;


  stn.adcp_sigwavehgt.date = dts;
  stn.adcp_sigwavehgt.data = [wvs.Hs]';
  stn.adcp_maxwavehgt.date = dts;
  stn.adcp_maxwavehgt.data = [wvs.Hmax]';
  stn.adcp_peakwaveper.date = dts;
  stn.adcp_peakwaveper.data = [wvs.Tp]';
  stn.adcp_peakwavedir.date = dts;
  stn.adcp_peakwavedir.data = [wvs.Dp]';
  stn.adcp_meanwaveper.date = dts;
  stn.adcp_meanwaveper.data = [wvs.Tmean]';

return;
