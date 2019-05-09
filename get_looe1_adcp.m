function stn = get_looe1_adcp(stn,doPlot)
%function stn = get_looe1_adcp([stn],[doPlot])
%
% Load multiple MAT files previously saved using WinADCP or similar, for ADCP
% (and bottom thermistor) data at AOML South Florida Program's Looe Key Reef
% spar. Do basic ADCP QC (percent good, sea-surface artifacts), and calculate
% barotropic averages, baroclinic residuals, and lowpass filter time series.
%
% See file THESIS/data/SFP/SERIAL_NUMBERS.txt for dates and deployment #s.
%
% Last Saved Time-stamp: <Mon 2012-04-09 12:07:34  Lew.Gramer>

  set_more off
  %DEBUG:  tic,

  %looepath = 'C:\Documents and Settings\gramer\My Documents\RSMAS\Coastal\thesis\data\SFP';
  datapath = get_thesis_path('../data');
  looepath = fullfile(datapath,'SFP');

  if ( ~exist('stn','var') )
    stn = [];
  end;
  if ( ~exist('doPlot','var') || isempty(doPlot) )
    doPlot = false;
  end;

  matfname = fullfile(datapath,'looe1_adcp.mat');
  if ( exist(matfname,'file') )
    disp(['Reloading pre-saved file ' matfname]);
    load(matfname,'result');

  else

    disp(['Processing raw ADCP MAT data...']);

    result.station_name = 'LOOE1';
    [result.lon,result.lat,result.depth] = get_station_coords(result.station_name);

    result.adcp_bin_heights = [];
    result.adcp_speed = struct('date',[],'data',[],'rawprof',[],'prof',[]);
    result.adcp_dir = struct('date',[],'data',[],'rawprof',[],'prof',[]);
    result.adcp_x = struct('date',[],'data',[],'rawprof',[],'prof',[]);
    result.adcp_l = struct('date',[],'data',[],'rawprof',[],'prof',[]);
    result.adcp_u = struct('date',[],'data',[],'rawprof',[],'prof',[]);
    result.adcp_v = struct('date',[],'data',[],'rawprof',[],'prof',[]);
    result.adcp_w = struct('date',[],'data',[],'rawprof',[],'prof',[]);
    result.adcp_err = struct('date',[],'data',[],'rawprof',[],'prof',[]);

    result.adcp_avg_eacnt = [];
    result.adcp_pct_good = [];

    adcp_eacnt = [];
    %adcp_pct_good = [];


    % 'LR_ADCP_DEP05.asc', ...
    % Looe_Key_ADCP_19OCT2007_thru_11JUN2008.mat
    % Looe_Key_ADCP_21JUN2008_thru_25JAN2009.mat
    % Looe_Ley_ADCP_25FEB2007_thru_19OCT2007.mat
    % LooeKey_ADCP_418_MAY2010.mat
    fnames = { ...
        'Looe_Ley_ADCP_05NOV2004_thru_26JAN2005.mat', ...
        'Looe_Key_ADCP_25MAR2005_thru_15SEP2005.mat', ...
        'Looe_Key_ADCP_15SEP2005_thru_16MAR2006.mat', ...
        'Looe_Key_ADCP_17MAR2006_thru_18JUL2006.mat', ...
        'Looe_Key_ADCP_18JUL2006_thru_24FEB2007.mat', ...
        'Looe_Ley_ADCP_25FEB2007_thru_19OCT2007.mat', ...
        'Looe_Key_ADCP_19OCT2007_thru_11JUN2008.mat', ...
        'Looe_Key_ADCP_21JUN2008_thru_25JAN2009.mat', ...
        'Looe_Key_ADCP_27MAY2009_thru_10FEB2010.mat', ...
             };

    for ix = 1:length(fnames)

      fname = fullfile(looepath,fnames{ix});
      %DEBUG:    disp(fname);
      looe = load(fname);

      ndts = length(looe.SerEnsembles);
      nbins = length(looe.SerBins);

      yrs = looe.SerYear;
      yrs(yrs < 70) = yrs(yrs < 70) + 2000;
      yrs(yrs < 100) = yrs(yrs < 100) + 1900;

      % HACKOLA! For some reason, the SerHour (and SerDay) in one or more of
      % these MAT files do not increase monotonically within days, but instead
      % have jumps up and back: and yet, NDTS == ROUNDN(DTS(end)-DTS(1))*24!
      %dts = datenum(yrs,looe.SerMon,looe.SerDay,looe.SerHour,looe.SerMin,looe.SerSec);
      % WORKAROUND: Assume no gaps, and that the first time stamp is correct.
      begdt = datenum(yrs(1),looe.SerMon(1),looe.SerDay(1),...
                      looe.SerHour(1),looe.SerMin(1),looe.SerSec(1));
      dts = begdt + [0:ndts-1]'*(1/24);

      result.adcp_bin_heights = looe.RDIBin1Mid + looe.RDIBinSize*(looe.SerBins-1);
      % result.adcp_bin_heights = -( ceil(max(result.adcp_bin_heights)) - result.adcp_bin_heights );

      endix = length(result.adcp_speed.date);

      result.adcp_seatemp.date(endix+1:endix+ndts,1) = dts;
      result.adcp_seatemp.data(endix+1:endix+ndts,1) = looe.AnT100thDeg ./ 100;

      result.adcp_speed.date(endix+1:endix+ndts,1) = dts;
      result.adcp_speed.prof(endix+1:endix+ndts,1:nbins) = looe.SerMagmmpersec ./ 1000;
      result.adcp_dir.date(endix+1:endix+ndts,1) = dts;
      result.adcp_dir.prof(endix+1:endix+ndts,1:nbins) = looe.SerDir10thDeg ./ 10;

      % Rotate U (V) current component to local cross (along) shore orientation
      u = looe.SerEmmpersec ./ 1000; % mm/s to m/s
      v = looe.SerNmmpersec ./ 1000;

      % Orientation of local isobath in degrees True
      % ori = 70; % Gradient-perpendicular direction from NGDC 3-arcsec bathymetry 
      % ori = 77; % Google Earth orientation estimate along a 1km leg
      ori = 73; % Published estimate in Lee & Williams (1999)

      result.adcp_isobath_orientation = ori;
      xshore = (cosd(ori)*u) - (sind(ori)*v);  % Dyslexics of the world untie
      lshore = (sind(ori)*u) + (cosd(ori)*v);

      result.adcp_x.date(endix+1:endix+ndts,1) = dts;
      result.adcp_x.prof(endix+1:endix+ndts,1:nbins) = xshore;
      result.adcp_l.date(endix+1:endix+ndts,1) = dts;
      result.adcp_l.prof(endix+1:endix+ndts,1:nbins) = lshore;
      result.adcp_u.date(endix+1:endix+ndts,1) = dts;
      result.adcp_u.prof(endix+1:endix+ndts,1:nbins) = u;
      result.adcp_v.date(endix+1:endix+ndts,1) = dts;
      result.adcp_v.prof(endix+1:endix+ndts,1:nbins) = v;
      result.adcp_w.date(endix+1:endix+ndts,1) = dts;
      result.adcp_w.prof(endix+1:endix+ndts,1:nbins) = looe.SerVmmpersec ./ 1000;
      result.adcp_err.date(endix+1:endix+ndts,1) = dts;
      result.adcp_err.prof(endix+1:endix+ndts,1:nbins) = looe.SerErmmpersec ./ 1000;

      result.adcp_pct_good.date(endix+1:endix+ndts) = dts;
      result.adcp_pct_good.prof(endix+1:endix+ndts,1:nbins) = looe.SerPG4;
      %result.adcp_ensemble_pct_good.date(endix+1:endix+ndts) = dts;
      %result.adcp_ensemble_pct_good.prof(endix+1:endix+ndts,1:nbins) = looe.SerPG3;

      adcp_eacnt(1,endix+1:endix+ndts,1:nbins) = looe.SerEA1cnt;
      adcp_eacnt(2,endix+1:endix+ndts,1:nbins) = looe.SerEA2cnt;
      adcp_eacnt(3,endix+1:endix+ndts,1:nbins) = looe.SerEA3cnt;
      % Workaround for weirdness in the 2007-2008 MAT file - ??
      if ( isfield(looe,'SerEA4cnt') )
        adcp_eacnt(4,endix+1:endix+ndts,1:nbins) = looe.SerEA4cnt;
      elseif ( isfield(looe,'SerEAAcnt') )
        adcp_eacnt(4,endix+1:endix+ndts,1:nbins) = looe.SerEAAcnt;
      end;

      looe = []; clear looe;

    end; %for ix = 1:length(fnames)

    % Save raw profiles for future reference, before surface-finding and QA/QC
    result.adcp_speed.rawprof = result.adcp_speed.prof;
    result.adcp_dir.rawprof = result.adcp_dir.prof;
    result.adcp_x.rawprof = result.adcp_x.prof;
    result.adcp_l.rawprof = result.adcp_l.prof;
    result.adcp_u.rawprof = result.adcp_u.prof;
    result.adcp_v.rawprof = result.adcp_v.prof;
    result.adcp_w.rawprof = result.adcp_w.prof;
    result.adcp_err.rawprof = result.adcp_err.prof;


    % Find surface using echo amplitudes with an algorithm from Ryan.Smith@noaa.gov
    % %average all four echo amplitudes over the length of the time series...
    % avg_EA = reshape((mean([SerEA1cnt(:)'; SerEA2cnt(:)'; SerEA3cnt(:)'; SerEA4cnt(:)'])),(length(SerEA1cnt)),36);
    % % 1) find the surface and calculate good bins...
    % % 2) look between bins 20 and 36 for the surface spike...
    % surface_bin=find(avg_EA(i,20:36)==max(avg_EA(i,20:36))) + 19; % add 19 to get true bin number...
    % % 3) beam angles are 20 deg, so cut last 6% plus 1 bin of profile to remove side lobe contam...
    % last_good_bin = surface_bin - (floor(0.06*surface_bin) + 1);
    %set good data to be all data with "percent good" higher than pgood;

    %DEBUG:  disp('tic'); tic,
    result.adcp_avg_eacnt.date = dts;
    result.adcp_avg_eacnt.prof = squeeze(mean(adcp_eacnt,1));
    max_eacnt = max(result.adcp_avg_eacnt.prof(:,20:36),[],2);

    % Loop over each timestamp, to find the bin number just below the surface spike
    % If I were a bigger MATLAB geek, I'd probably figure out how to remove this loop!
    for ix=1:size(result.adcp_avg_eacnt.prof,1)
      %DEBUG:    if (mod(ix,1000)==1); disp(datestr(result.adcp_u.date(ix))); end;
      result.adcp_sfc_bin(ix,1) = find(result.adcp_avg_eacnt.prof(ix,20:36)==max_eacnt(ix),1) + 19;
      result.adcp_sfc_height(ix,1) = result.adcp_bin_heights(result.adcp_sfc_bin(ix));

      % Beam angles 20o: side lobe contamination in last 6% + 1 bin of profile
      last_good_bin = result.adcp_sfc_bin(ix) - (floor(0.06*result.adcp_sfc_bin(ix)) + 1);
      result.adcp_last_good_bin(ix,1) = last_good_bin;
      result.adcp_speed.prof(ix,last_good_bin+1:end) = nan;
      result.adcp_dir.prof(ix,last_good_bin+1:end) = nan;
      result.adcp_x.prof(ix,last_good_bin+1:end) = nan;
      result.adcp_l.prof(ix,last_good_bin+1:end) = nan;
      result.adcp_u.prof(ix,last_good_bin+1:end) = nan;
      result.adcp_v.prof(ix,last_good_bin+1:end) = nan;
      result.adcp_w.prof(ix,last_good_bin+1:end) = nan;
      result.adcp_err.prof(ix,last_good_bin+1:end) = nan;
    end;
    %DEBUG:  toc,

    %DEBUG:
    disp('QC');
    % Replace supposedly spurious returns with NaNs
    badix = find(result.adcp_pct_good.prof(:) < 50);
    result.adcp_speed.prof(badix) = nan;
    result.adcp_dir.prof(badix) = nan;
    result.adcp_x.prof(badix) = nan;
    result.adcp_l.prof(badix) = nan;
    result.adcp_u.prof(badix) = nan;
    result.adcp_v.prof(badix) = nan;
    result.adcp_w.prof(badix) = nan;
    result.adcp_err.prof(badix) = nan;

    % After removal of spurious returns, also do barotropic averages
    result.adcp_speed.data = nanmean(result.adcp_speed.prof,2);
    result.adcp_x.data = nanmean(result.adcp_x.prof,2);
    result.adcp_l.data = nanmean(result.adcp_l.prof,2);
    result.adcp_u.data = nanmean(result.adcp_u.prof,2);
    result.adcp_v.data = nanmean(result.adcp_v.prof,2);
    result.adcp_w.data = nanmean(result.adcp_w.prof,2);
    result.adcp_err.data = nanmean(result.adcp_err.prof,2);
    % Don't forget to use vectorial averaging
    result.adcp_dir.data = uv_to_dir_curr(result.adcp_u.data,result.adcp_v.data);

    % Calculate 'baroclinic' currents by removing barotropic averages
    result.adcp_baroclinic_speed.date = result.adcp_speed.date;
    result.adcp_baroclinic_speed.prof = result.adcp_speed.prof - repmat(result.adcp_speed.data,[1 nbins]);
    result.adcp_baroclinic_x.date = result.adcp_x.date;
    result.adcp_baroclinic_x.prof = result.adcp_x.prof - repmat(result.adcp_x.data,[1 nbins]);
    result.adcp_baroclinic_l.date = result.adcp_l.date;
    result.adcp_baroclinic_l.prof = result.adcp_l.prof - repmat(result.adcp_l.data,[1 nbins]);
    result.adcp_baroclinic_u.date = result.adcp_u.date;
    result.adcp_baroclinic_u.prof = result.adcp_u.prof - repmat(result.adcp_u.data,[1 nbins]);
    result.adcp_baroclinic_v.date = result.adcp_v.date;
    result.adcp_baroclinic_v.prof = result.adcp_v.prof - repmat(result.adcp_v.data,[1 nbins]);
    result.adcp_baroclinic_dir.date = result.adcp_dir.date;
    result.adcp_baroclinic_dir.prof = uv_to_dir_curr(result.adcp_baroclinic_u.prof,result.adcp_baroclinic_v.prof);


    %DEBUG:
    disp('3HLP/40HLP');
    % Also low-pass filter using LANCZOSFILTER (v.) and GAP_FILTER_TS (v.)

    [ig,result.adcp_speed_3hlp.data] = lanczos_ts(result.adcp_speed.date,result.adcp_speed.data);
    [ig,result.adcp_dir_3hlp.data] = lanczos_ts(result.adcp_dir.date,result.adcp_dir.data);
    [ig,result.adcp_x_3hlp.data] = lanczos_ts(result.adcp_x.date,result.adcp_x.data);
    [ig,result.adcp_l_3hlp.data] = lanczos_ts(result.adcp_l.date,result.adcp_l.data);
    [ig,result.adcp_u_3hlp.data] = lanczos_ts(result.adcp_u.date,result.adcp_u.data);
    [ig,result.adcp_v_3hlp.data] = lanczos_ts(result.adcp_v.date,result.adcp_v.data);
    [ig,result.adcp_w_3hlp.data] = lanczos_ts(result.adcp_w.date,result.adcp_w.data);
    [ig,result.adcp_err_3hlp.data] = lanczos_ts(result.adcp_err.date,result.adcp_err.data);

    [ig,result.adcp_speed_40hlp.data] = gap_filter_ts(result.adcp_speed.date,result.adcp_speed.data,40);
    [ig,result.adcp_dir_40hlp.data] = gap_filter_ts(result.adcp_dir.date,result.adcp_dir.data,40);
    [ig,result.adcp_x_40hlp.data] = gap_filter_ts(result.adcp_x.date,result.adcp_x.data,40);
    [ig,result.adcp_l_40hlp.data] = gap_filter_ts(result.adcp_l.date,result.adcp_l.data,40);
    [ig,result.adcp_u_40hlp.data] = gap_filter_ts(result.adcp_u.date,result.adcp_u.data,40);
    [ig,result.adcp_v_40hlp.data] = gap_filter_ts(result.adcp_v.date,result.adcp_v.data,40);
    [ig,result.adcp_w_40hlp.data] = gap_filter_ts(result.adcp_w.date,result.adcp_w.data,40);
    [ig,result.adcp_err_40hlp.data] = gap_filter_ts(result.adcp_err.date,result.adcp_err.data,40);

    dts = [];
    % Again, this code is ripe for some kind of genius vectorialization
    for ix=1:nbins
      %DEBUG:
      disp(ix);
      % Shortcut for bins with no data (e.g., "above-surface" bins)
      if ( all(isnan(result.adcp_dir.prof(:,ix))) && ~isempty(dts) )
        result.adcp_speed_3hlp.prof(:,ix) = nan;
        result.adcp_dir_3hlp.prof(:,ix) = nan;
        result.adcp_x_3hlp.prof(:,ix) = nan;
        result.adcp_l_3hlp.prof(:,ix) = nan;
        result.adcp_u_3hlp.prof(:,ix) = nan;
        result.adcp_v_3hlp.prof(:,ix) = nan;
        result.adcp_w_3hlp.prof(:,ix) = nan;
        result.adcp_err_3hlp.prof(:,ix) = nan;

        result.adcp_baroclinic_speed_3hlp.prof(:,ix) = nan;
        result.adcp_baroclinic_dir_3hlp.prof(:,ix) = nan;
        result.adcp_baroclinic_x_3hlp.prof(:,ix) = nan;
        result.adcp_baroclinic_l_3hlp.prof(:,ix) = nan;
        result.adcp_baroclinic_u_3hlp.prof(:,ix) = nan;
        result.adcp_baroclinic_v_3hlp.prof(:,ix) = nan;

        result.adcp_speed_40hlp.prof(:,ix) = nan;
        result.adcp_dir_40hlp.prof(:,ix) = nan;
        result.adcp_x_40hlp.prof(:,ix) = nan;
        result.adcp_l_40hlp.prof(:,ix) = nan;
        result.adcp_u_40hlp.prof(:,ix) = nan;
        result.adcp_v_40hlp.prof(:,ix) = nan;
        result.adcp_w_40hlp.prof(:,ix) = nan;
        result.adcp_err_40hlp.prof(:,ix) = nan;

        result.adcp_baroclinic_speed_40hlp.prof(:,ix) = nan;
        result.adcp_baroclinic_dir_40hlp.prof(:,ix) = nan;
        result.adcp_baroclinic_x_40hlp.prof(:,ix) = nan;
        result.adcp_baroclinic_l_40hlp.prof(:,ix) = nan;
        result.adcp_baroclinic_u_40hlp.prof(:,ix) = nan;
        result.adcp_baroclinic_v_40hlp.prof(:,ix) = nan;

      else

        [dts,dat] = lanczos_ts(result.adcp_speed.date,result.adcp_speed.prof(:,ix));
        ndts = length(dts);

        result.adcp_speed_3hlp.date = dts;
        result.adcp_speed_3hlp.prof(1:ndts,ix) = dat;
        [dts,dat] = lanczos_ts(result.adcp_dir.date,result.adcp_dir.prof(:,ix));
        result.adcp_dir_3hlp.date = dts;
        result.adcp_dir_3hlp.prof(1:ndts,ix) = dat;
        [dts,dat] = lanczos_ts(result.adcp_x.date,result.adcp_x.prof(:,ix));
        result.adcp_x_3hlp.date = dts;
        result.adcp_x_3hlp.prof(1:ndts,ix) = dat;
        [dts,dat] = lanczos_ts(result.adcp_l.date,result.adcp_l.prof(:,ix));
        result.adcp_l_3hlp.date = dts;
        result.adcp_l_3hlp.prof(1:ndts,ix) = dat;
        [dts,dat] = lanczos_ts(result.adcp_u.date,result.adcp_u.prof(:,ix));
        result.adcp_u_3hlp.date = dts;
        result.adcp_u_3hlp.prof(1:ndts,ix) = dat;
        [dts,dat] = lanczos_ts(result.adcp_v.date,result.adcp_v.prof(:,ix));
        result.adcp_v_3hlp.date = dts;
        result.adcp_v_3hlp.prof(1:ndts,ix) = dat;
        [dts,dat] = lanczos_ts(result.adcp_w.date,result.adcp_w.prof(:,ix));
        result.adcp_w_3hlp.date = dts;
        result.adcp_w_3hlp.prof(1:ndts,ix) = dat;
        [dts,dat] = lanczos_ts(result.adcp_err.date,result.adcp_err.prof(:,ix));
        result.adcp_err_3hlp.date = dts;
        result.adcp_err_3hlp.prof(1:ndts,ix) = dat;

        [dts,dat] = lanczos_ts(result.adcp_baroclinic_speed.date,result.adcp_baroclinic_speed.prof(:,ix));
        result.adcp_baroclinic_speed_3hlp.date = dts;
        result.adcp_baroclinic_speed_3hlp.prof(1:ndts,ix) = dat;
        [dts,dat] = lanczos_ts(result.adcp_baroclinic_dir.date,result.adcp_baroclinic_dir.prof(:,ix));
        result.adcp_baroclinic_dir_3hlp.date = dts;
        result.adcp_baroclinic_dir_3hlp.prof(1:ndts,ix) = dat;
        [dts,dat] = lanczos_ts(result.adcp_baroclinic_x.date,result.adcp_baroclinic_x.prof(:,ix));
        result.adcp_baroclinic_x_3hlp.date = dts;
        result.adcp_baroclinic_x_3hlp.prof(1:ndts,ix) = dat;
        [dts,dat] = lanczos_ts(result.adcp_baroclinic_l.date,result.adcp_baroclinic_l.prof(:,ix));
        result.adcp_baroclinic_l_3hlp.date = dts;
        result.adcp_baroclinic_l_3hlp.prof(1:ndts,ix) = dat;
        [dts,dat] = lanczos_ts(result.adcp_baroclinic_u.date,result.adcp_baroclinic_u.prof(:,ix));
        result.adcp_baroclinic_u_3hlp.date = dts;
        result.adcp_baroclinic_u_3hlp.prof(1:ndts,ix) = dat;
        [dts,dat] = lanczos_ts(result.adcp_baroclinic_v.date,result.adcp_baroclinic_v.prof(:,ix));
        result.adcp_baroclinic_v_3hlp.date = dts;
        result.adcp_baroclinic_v_3hlp.prof(1:ndts,ix) = dat;

        [dts40,dat40] = gap_filter_ts(result.adcp_speed.date,result.adcp_speed.prof(:,ix),40);
        ndts40 = length(dts40);

        result.adcp_speed_40hlp.date = dts40;
        result.adcp_speed_40hlp.prof(1:ndts40,ix) = dat40;
        [dts40,dat40] = gap_filter_ts(result.adcp_dir.date,result.adcp_dir.prof(:,ix),40);
        result.adcp_dir_40hlp.date = dts40;
        result.adcp_dir_40hlp.prof(1:ndts40,ix) = dat40;
        [dts40,dat40] = gap_filter_ts(result.adcp_x.date,result.adcp_x.prof(:,ix),40);
        result.adcp_x_40hlp.date = dts40;
        result.adcp_x_40hlp.prof(1:ndts40,ix) = dat40;
        [dts40,dat40] = gap_filter_ts(result.adcp_l.date,result.adcp_l.prof(:,ix),40);
        result.adcp_l_40hlp.date = dts40;
        result.adcp_l_40hlp.prof(1:ndts40,ix) = dat40;
        [dts40,dat40] = gap_filter_ts(result.adcp_u.date,result.adcp_u.prof(:,ix),40);
        result.adcp_u_40hlp.date = dts40;
        result.adcp_u_40hlp.prof(1:ndts40,ix) = dat40;
        [dts40,dat40] = gap_filter_ts(result.adcp_v.date,result.adcp_v.prof(:,ix),40);
        result.adcp_v_40hlp.date = dts40;
        result.adcp_v_40hlp.prof(1:ndts40,ix) = dat40;
        [dts40,dat40] = gap_filter_ts(result.adcp_w.date,result.adcp_w.prof(:,ix),40);
        result.adcp_w_40hlp.date = dts40;
        result.adcp_w_40hlp.prof(1:ndts40,ix) = dat40;
        [dts40,dat40] = gap_filter_ts(result.adcp_err.date,result.adcp_err.prof(:,ix),40);
        result.adcp_err_40hlp.date = dts40;
        result.adcp_err_40hlp.prof(1:ndts40,ix) = dat40;

        [dts40,dat40] = gap_filter_ts(result.adcp_baroclinic_speed.date,result.adcp_baroclinic_speed.prof(:,ix),40);
        result.adcp_baroclinic_speed_40hlp.date = dts40;
        result.adcp_baroclinic_speed_40hlp.prof(1:ndts40,ix) = dat40;
        [dts40,dat40] = gap_filter_ts(result.adcp_baroclinic_dir.date,result.adcp_baroclinic_dir.prof(:,ix),40);
        result.adcp_baroclinic_dir_40hlp.date = dts40;
        result.adcp_baroclinic_dir_40hlp.prof(1:ndts40,ix) = dat40;
        [dts40,dat40] = gap_filter_ts(result.adcp_baroclinic_x.date,result.adcp_baroclinic_x.prof(:,ix),40);
        result.adcp_baroclinic_x_40hlp.date = dts40;
        result.adcp_baroclinic_x_40hlp.prof(1:ndts40,ix) = dat40;
        [dts40,dat40] = gap_filter_ts(result.adcp_baroclinic_l.date,result.adcp_baroclinic_l.prof(:,ix),40);
        result.adcp_baroclinic_l_40hlp.date = dts40;
        result.adcp_baroclinic_l_40hlp.prof(1:ndts40,ix) = dat40;
        [dts40,dat40] = gap_filter_ts(result.adcp_baroclinic_u.date,result.adcp_baroclinic_u.prof(:,ix),40);
        result.adcp_baroclinic_u_40hlp.date = dts40;
        result.adcp_baroclinic_u_40hlp.prof(1:ndts40,ix) = dat40;
        [dts40,dat40] = gap_filter_ts(result.adcp_baroclinic_v.date,result.adcp_baroclinic_v.prof(:,ix),40);
        result.adcp_baroclinic_v_40hlp.date = dts40;
        result.adcp_baroclinic_v_40hlp.prof(1:ndts40,ix) = dat40;

      end; %if ( all(isnan(result.adcp_dir.prof(:,ix))) && ~isempty(dts) ) else

    end; %for ix=1:nbins

    disp(['Saving result to file ' matfname]);
    save(matfname,'result');

  end; % if exist(matfname) else

  flds = fieldnames(result);
  for ix = 1:length(flds)
    fld = flds{ix};
    stn.(fld) = result.(fld);
  end;
  result = []; clear result;



  %% Do multi-bin averages: "surface", "mid", "bottom", and 3 bins around MicroCAT


  % Near-surface (top 6 bins, 4.5 m below any side-lobe contamination)
  stn.adcp_sfc_x.date = stn.adcp_x.date;
  stn.adcp_sfc_l.date = stn.adcp_l.date;
  stn.adcp_sfc_u.date = stn.adcp_u.date;
  stn.adcp_sfc_v.date = stn.adcp_v.date;
  stn.adcp_sfc_w.date = stn.adcp_w.date;
  stn.adcp_baroclinic_sfc_x.date = stn.adcp_baroclinic_x.date;
  stn.adcp_baroclinic_sfc_l.date = stn.adcp_baroclinic_l.date;
  stn.adcp_baroclinic_sfc_u.date = stn.adcp_baroclinic_u.date;
  stn.adcp_baroclinic_sfc_v.date = stn.adcp_baroclinic_v.date;
  for ix = 1:numel(stn.adcp_sfc_u.date)
    tix = stn.adcp_last_good_bin(ix);
    sfcix = tix-5:tix;
    stn.adcp_sfc_x.data(ix,1) = nanmean(stn.adcp_x.prof(ix,sfcix));
    stn.adcp_sfc_l.data(ix,1) = nanmean(stn.adcp_l.prof(ix,sfcix));
    stn.adcp_sfc_u.data(ix,1) = nanmean(stn.adcp_u.prof(ix,sfcix));
    stn.adcp_sfc_v.data(ix,1) = nanmean(stn.adcp_v.prof(ix,sfcix));
    stn.adcp_sfc_w.data(ix,1) = nanmean(stn.adcp_w.prof(ix,sfcix));
    stn.adcp_baroclinic_sfc_x.data(ix,1) = nanmean(stn.adcp_baroclinic_x.prof(ix,sfcix));
    stn.adcp_baroclinic_sfc_l.data(ix,1) = nanmean(stn.adcp_baroclinic_l.prof(ix,sfcix));
    stn.adcp_baroclinic_sfc_u.data(ix,1) = nanmean(stn.adcp_baroclinic_u.prof(ix,sfcix));
    stn.adcp_baroclinic_sfc_v.data(ix,1) = nanmean(stn.adcp_baroclinic_v.prof(ix,sfcix));
  end;
  stn.adcp_sfc_speed.date = stn.adcp_sfc_u.date;
  stn.adcp_sfc_speed.data = uv_to_spd(stn.adcp_sfc_u.data,stn.adcp_sfc_v.data);
  stn.adcp_sfc_dir.date = stn.adcp_sfc_u.date;
  stn.adcp_sfc_dir.data = uv_to_dir_curr(stn.adcp_sfc_u.data,stn.adcp_sfc_v.data);

  stn.adcp_baroclinic_sfc_speed.date = stn.adcp_baroclinic_sfc_u.date;
  stn.adcp_baroclinic_sfc_speed.data = uv_to_spd(stn.adcp_baroclinic_sfc_u.data,stn.adcp_baroclinic_sfc_v.data);
  stn.adcp_baroclinic_sfc_dir.date = stn.adcp_baroclinic_sfc_u.date;
  stn.adcp_baroclinic_sfc_dir.data = uv_to_dir_curr(stn.adcp_baroclinic_sfc_u.data,stn.adcp_baroclinic_sfc_v.data);


  % Mid-column (bins 11-16, from 9.3 to 13.8m above bottom
  midix = 11:16;
  stn.adcp_mid_x.date = stn.adcp_x.date;
  stn.adcp_mid_x.data = nanmean(stn.adcp_x.prof(:,midix),2);
  stn.adcp_mid_l.date = stn.adcp_l.date;
  stn.adcp_mid_l.data = nanmean(stn.adcp_l.prof(:,midix),2);
  stn.adcp_mid_u.date = stn.adcp_u.date;
  stn.adcp_mid_u.data = nanmean(stn.adcp_u.prof(:,midix),2);
  stn.adcp_mid_v.date = stn.adcp_v.date;
  stn.adcp_mid_v.data = nanmean(stn.adcp_v.prof(:,midix),2);
  stn.adcp_mid_w.date = stn.adcp_w.date;
  stn.adcp_mid_w.data = nanmean(stn.adcp_w.prof(:,midix),2);

  stn.adcp_mid_speed.date = stn.adcp_mid_u.date;
  stn.adcp_mid_speed.data = uv_to_spd(stn.adcp_mid_u.data,stn.adcp_mid_v.data);
  stn.adcp_mid_dir.date = stn.adcp_mid_u.date;
  stn.adcp_mid_dir.data = uv_to_dir_curr(stn.adcp_mid_u.data,stn.adcp_mid_v.data);

  stn.adcp_baroclinic_mid_x.date = stn.adcp_baroclinic_x.date;
  stn.adcp_baroclinic_mid_x.data = nanmean(stn.adcp_baroclinic_x.prof(:,midix),2);
  stn.adcp_baroclinic_mid_l.date = stn.adcp_baroclinic_l.date;
  stn.adcp_baroclinic_mid_l.data = nanmean(stn.adcp_baroclinic_l.prof(:,midix),2);
  stn.adcp_baroclinic_mid_u.date = stn.adcp_baroclinic_u.date;
  stn.adcp_baroclinic_mid_u.data = nanmean(stn.adcp_baroclinic_u.prof(:,midix),2);
  stn.adcp_baroclinic_mid_v.date = stn.adcp_baroclinic_v.date;
  stn.adcp_baroclinic_mid_v.data = nanmean(stn.adcp_baroclinic_v.prof(:,midix),2);

  stn.adcp_baroclinic_mid_speed.date = stn.adcp_baroclinic_mid_u.date;
  stn.adcp_baroclinic_mid_speed.data = uv_to_spd(stn.adcp_baroclinic_mid_u.data,stn.adcp_baroclinic_mid_v.data);
  stn.adcp_baroclinic_mid_dir.date = stn.adcp_baroclinic_mid_u.date;
  stn.adcp_baroclinic_mid_dir.data = uv_to_dir_curr(stn.adcp_baroclinic_mid_u.data,stn.adcp_baroclinic_mid_v.data);


  % Three bins surrounding MicroCAT depth (5m below mean surface, ~bin 21)
  mcdix = 20:22;
  badix = find(stn.adcp_last_good_bin<max(mcdix));
  stn.adcp_mcd_x.date = stn.adcp_x.date;
  stn.adcp_mcd_x.data = nanmean(stn.adcp_x.prof(:,mcdix),2);
  stn.adcp_mcd_x.data(badix) = nan;
  stn.adcp_mcd_l.date = stn.adcp_l.date;
  stn.adcp_mcd_l.data = nanmean(stn.adcp_l.prof(:,mcdix),2);
  stn.adcp_mcd_l.data(badix) = nan;
  stn.adcp_mcd_u.date = stn.adcp_u.date;
  stn.adcp_mcd_u.data = nanmean(stn.adcp_u.prof(:,mcdix),2);
  stn.adcp_mcd_u.data(badix) = nan;
  stn.adcp_mcd_v.date = stn.adcp_v.date;
  stn.adcp_mcd_v.data = nanmean(stn.adcp_v.prof(:,mcdix),2);
  stn.adcp_mcd_v.data(badix) = nan;
  stn.adcp_mcd_w.date = stn.adcp_w.date;
  stn.adcp_mcd_w.data = nanmean(stn.adcp_w.prof(:,mcdix),2);
  stn.adcp_mcd_w.data(badix) = nan;

  stn.adcp_mcd_speed.date = stn.adcp_mcd_u.date;
  stn.adcp_mcd_speed.data = uv_to_spd(stn.adcp_mcd_u.data,stn.adcp_mcd_v.data);
  stn.adcp_mcd_dir.date = stn.adcp_mcd_u.date;
  stn.adcp_mcd_dir.data = uv_to_dir_curr(stn.adcp_mcd_u.data,stn.adcp_mcd_v.data);

  stn.adcp_baroclinic_mcd_x.date = stn.adcp_baroclinic_x.date;
  stn.adcp_baroclinic_mcd_x.data = nanmean(stn.adcp_baroclinic_x.prof(:,mcdix),2);
  stn.adcp_baroclinic_mcd_x.data(badix) = nan;
  stn.adcp_baroclinic_mcd_l.date = stn.adcp_baroclinic_l.date;
  stn.adcp_baroclinic_mcd_l.data = nanmean(stn.adcp_baroclinic_l.prof(:,mcdix),2);
  stn.adcp_baroclinic_mcd_l.data(badix) = nan;
  stn.adcp_baroclinic_mcd_u.date = stn.adcp_baroclinic_u.date;
  stn.adcp_baroclinic_mcd_u.data = nanmean(stn.adcp_baroclinic_u.prof(:,mcdix),2);
  stn.adcp_baroclinic_mcd_u.data(badix) = nan;
  stn.adcp_baroclinic_mcd_v.date = stn.adcp_baroclinic_v.date;
  stn.adcp_baroclinic_mcd_v.data = nanmean(stn.adcp_baroclinic_v.prof(:,mcdix),2);
  stn.adcp_baroclinic_mcd_v.data(badix) = nan;

  stn.adcp_baroclinic_mcd_speed.date = stn.adcp_baroclinic_mcd_u.date;
  stn.adcp_baroclinic_mcd_speed.data = uv_to_spd(stn.adcp_baroclinic_mcd_u.data,stn.adcp_baroclinic_mcd_v.data);
  stn.adcp_baroclinic_mcd_dir.date = stn.adcp_baroclinic_mcd_u.date;
  stn.adcp_baroclinic_mcd_dir.data = uv_to_dir_curr(stn.adcp_baroclinic_mcd_u.data,stn.adcp_baroclinic_mcd_v.data);


  % Near-bottom (first 6 bins, 4.5 m - above blanking area)
  btmix = 1:6;
  stn.adcp_btm_x.date = stn.adcp_x.date;
  stn.adcp_btm_x.data = nanmean(stn.adcp_x.prof(:,btmix),2);
  stn.adcp_btm_l.date = stn.adcp_l.date;
  stn.adcp_btm_l.data = nanmean(stn.adcp_l.prof(:,btmix),2);
  stn.adcp_btm_u.date = stn.adcp_u.date;
  stn.adcp_btm_u.data = nanmean(stn.adcp_u.prof(:,btmix),2);
  stn.adcp_btm_v.date = stn.adcp_v.date;
  stn.adcp_btm_v.data = nanmean(stn.adcp_v.prof(:,btmix),2);
  stn.adcp_btm_w.date = stn.adcp_w.date;
  stn.adcp_btm_w.data = nanmean(stn.adcp_w.prof(:,btmix),2);

  stn.adcp_btm_speed.date = stn.adcp_btm_u.date;
  stn.adcp_btm_speed.data = uv_to_spd(stn.adcp_btm_u.data,stn.adcp_btm_v.data);
  stn.adcp_btm_dir.date = stn.adcp_btm_u.date;
  stn.adcp_btm_dir.data = uv_to_dir_curr(stn.adcp_btm_u.data,stn.adcp_btm_v.data);

  stn.adcp_baroclinic_btm_x.date = stn.adcp_baroclinic_x.date;
  stn.adcp_baroclinic_btm_x.data = nanmean(stn.adcp_baroclinic_x.prof(:,btmix),2);
  stn.adcp_baroclinic_btm_l.date = stn.adcp_baroclinic_l.date;
  stn.adcp_baroclinic_btm_l.data = nanmean(stn.adcp_baroclinic_l.prof(:,btmix),2);
  stn.adcp_baroclinic_btm_u.date = stn.adcp_baroclinic_u.date;
  stn.adcp_baroclinic_btm_u.data = nanmean(stn.adcp_baroclinic_u.prof(:,btmix),2);
  stn.adcp_baroclinic_btm_v.date = stn.adcp_baroclinic_v.date;
  stn.adcp_baroclinic_btm_v.data = nanmean(stn.adcp_baroclinic_v.prof(:,btmix),2);

  stn.adcp_baroclinic_btm_speed.date = stn.adcp_baroclinic_btm_u.date;
  stn.adcp_baroclinic_btm_speed.data = uv_to_spd(stn.adcp_baroclinic_btm_u.data,stn.adcp_baroclinic_btm_v.data);
  stn.adcp_baroclinic_btm_dir.date = stn.adcp_baroclinic_btm_u.date;
  stn.adcp_baroclinic_btm_dir.data = uv_to_dir_curr(stn.adcp_baroclinic_btm_u.data,stn.adcp_baroclinic_btm_v.data);


  if ( doPlot )
    plot_adcp_bins(stn,25,2);
  end;

  %DEBUG:  toc,
  set_more;

return;
