function [klgf1,mrtf1] = an_ncore
%function [klgf1,mrtf1] = an_ncore
%
% Process currents/temp. data from the 2007-2008 NCORE experiment, provided
% to L Gramer courtesy of TN Lee and S Sponaugle, by K Shulzitski in 2011.
% Returns two STATION structs KLGF1 (Key Largo site offshore of French Reef)
% and MRTF1 (Marathon site offshore of Tennessee Reef) with 'shallow' (7 or 8
% m, resp.) and 'deep' (23 m) sensor time-series sub-structures.
%
% Last Saved Time-stamp: <Wed 2013-01-16 14:48:19 Eastern Standard Time gramer>

  matfname=fullfile(get_thesis_path('../data'),'NCORE.mat');

  if ( exist(matfname,'file') )

    disp(['Loading from ',matfname]);
    load(matfname,'stn');

  else

    fnames = {'2007CurrentandTemp_RawData.xls','2008CurrentandTemp_RawData.xlsx'};
    all_shtnms = {...
        {'k1lar101','k2lar101','k3mar101','k4mar101'},...
        {'KLT0801','KLB0801','MKT0801','MKB0801'},...
                 };

    for fnix = 1:numel(fnames)
      fname = fnames{fnix};
      shtnms = all_shtnms{fnix};
      x = importdata(fullfile(get_thesis_path('../data'),'NCORE',fname));

      for cshtnm=shtnms;
        shtnm = cshtnm{:};
        disp([fname,'.',shtnm]);
        sht=x.data.(shtnm);
        dts = datenum(sht(:,[3,1,2,4,5,6]));
        u = sht(:,9);
        v = sht(:,10);
        w = sht(:,11);
        P = sht(:,20);
        T = sht(:,21);
        sht=[]; clear sht;
        stn.(shtnm).cm_u_raw.date = dts;
        stn.(shtnm).cm_u_raw.data = u;
        stn.(shtnm).cm_v_raw.date = dts;
        stn.(shtnm).cm_v_raw.data = v;
        stn.(shtnm).cm_w_raw.date = dts;
        stn.(shtnm).cm_w_raw.data = w;
        stn.(shtnm).cm_seapres_raw.date = dts;
        stn.(shtnm).cm_seapres_raw.data = P;
        stn.(shtnm).cm_seatemp_raw.date = dts;
        stn.(shtnm).cm_seatemp_raw.data = T;
      end; %for cshtnm=shtnms;
      x=[]; clear x;
    end; %for fnix = 1:numel(fnames)

    disp(['Saving to ',matfname]);
    save(matfname,'stn');

  end; %if ( exist(matfname,'file') ) else


  stn.k1lar101 = station_reorient_vectors(stn.k1lar101,55,'cm_u_raw','cm_v_raw');
  stn.k2lar101 = station_reorient_vectors(stn.k2lar101,55,'cm_u_raw','cm_v_raw');
  stn.k3mar101 = station_reorient_vectors(stn.k3mar101,55,'cm_u_raw','cm_v_raw');
  stn.k4mar101 = station_reorient_vectors(stn.k4mar101,55,'cm_u_raw','cm_v_raw');

  stn.KLT0801 = station_reorient_vectors(stn.KLT0801,55,'cm_u_raw','cm_v_raw');
  stn.KLB0801 = station_reorient_vectors(stn.KLB0801,55,'cm_u_raw','cm_v_raw');
  stn.MKT0801 = station_reorient_vectors(stn.MKT0801,55,'cm_u_raw','cm_v_raw');
  stn.MKB0801 = station_reorient_vectors(stn.MKB0801,55,'cm_u_raw','cm_v_raw');


  % QA/QC
  stn.k1lar101 = an_ncore_qa(stn.k1lar101);
  stn.k2lar101 = an_ncore_qa(stn.k2lar101);
  stn.k3mar101 = an_ncore_qa(stn.k3mar101);
  stn.k4mar101 = an_ncore_qa(stn.k4mar101);

  stn.KLT0801 = an_ncore_qa(stn.KLT0801);
  stn.KLB0801 = an_ncore_qa(stn.KLB0801);
  stn.MKT0801 = an_ncore_qa(stn.MKT0801);
  stn.MKB0801 = an_ncore_qa(stn.MKB0801);

  stn.k1lar101 = station_reorient_vectors(stn.k1lar101,55,'cm_u_qc','cm_v_qc');
  stn.k2lar101 = station_reorient_vectors(stn.k2lar101,55,'cm_u_qc','cm_v_qc');
  stn.k3mar101 = station_reorient_vectors(stn.k3mar101,55,'cm_u_qc','cm_v_qc');
  stn.k4mar101 = station_reorient_vectors(stn.k4mar101,55,'cm_u_qc','cm_v_qc');

  stn.KLT0801 = station_reorient_vectors(stn.KLT0801,55,'cm_u_qc','cm_v_qc');
  stn.KLB0801 = station_reorient_vectors(stn.KLB0801,55,'cm_u_qc','cm_v_qc');
  stn.MKT0801 = station_reorient_vectors(stn.MKT0801,55,'cm_u_qc','cm_v_qc');
  stn.MKB0801 = station_reorient_vectors(stn.MKB0801,55,'cm_u_qc','cm_v_qc');


  % On 1/4/2013 4:17 PM, Kathryn Shulzitski wrote:
  % > 25°01.883' 	080°20.881'
  % > Key Largo meter
  % > 
  % > Tennessee Reef mooring coordinates:
  % > 2444.418 	8046.582

  klgf1.station_name = 'klgf1';
  % % NOTE: French Reef coordinates
  % klgf1.lat =  25.034283;
  % klgf1.lon = -80.348217;

  % Per K. Shulzitski
  klgf1.lat =  25.0313833;
  klgf1.lon = -80.3480167;

  klgf1.long_desc = 'Key Largo Current and Temperature Data from 7 and 23 m (meters located just offshore of French Reef)';
  klgf1.shallow_depth = 7;
  klgf1.depth = 23;


  mrtf1.station_name = 'mrtf1';
  % % NOTE: Tennessee Reef coordinates
  % mrtf1.lat = 24.746053;
  % mrtf1.lon = -80.782375;

  % Per K. Shulzitski
  mrtf1.lat = 24.7403000;
  mrtf1.lon = -80.7763667;

  mrtf1.long_desc = 'Marathon Current and Temperature Data from 8 and 23 m (meters located just offshore of Tennessee Reef)';
  mrtf1.shallow_depth = 8;
  mrtf1.depth = 23;


  for cfld={'u','v','w','seapres','seatemp','xshore','lshore'};

    fld = cfld{:};
    rfld = ['cm_',fld,'_raw'];
    qfld = ['cm_',fld,'_qc'];

    sfld = ['cm_shallow_',fld];
    dfld = ['cm_deep_',fld];
    srfld = [sfld,'_raw'];
    drfld = [dfld,'_raw'];
    sqfld = [sfld,'_10min'];
    dqfld = [dfld,'_10min'];

    klgf1.(srfld) = stn.k1lar101.(rfld);
    n = numel(stn.KLT0801.(rfld).date);
    klgf1.(srfld).date(end+1:end+n,1) = stn.KLT0801.(rfld).date;
    klgf1.(srfld).data(end+1:end+n,1) = stn.KLT0801.(rfld).data;
    % QA'd 10-minute data
    klgf1.(sqfld) = stn.k1lar101.(qfld);
    n = numel(stn.KLT0801.(qfld).date);
    klgf1.(sqfld).date(end+1:end+n,1) = stn.KLT0801.(qfld).date;
    klgf1.(sqfld).data(end+1:end+n,1) = stn.KLT0801.(qfld).data;

    klgf1.(drfld) = stn.k2lar101.(rfld);
    n = numel(stn.KLB0801.(rfld).date);
    klgf1.(drfld).date(end+1:end+n,1) = stn.KLB0801.(rfld).date;
    klgf1.(drfld).data(end+1:end+n,1) = stn.KLB0801.(rfld).data;
    % QA'd 10-minute data
    klgf1.(dqfld) = stn.k2lar101.(qfld);
    n = numel(stn.KLB0801.(qfld).date);
    klgf1.(dqfld).date(end+1:end+n,1) = stn.KLB0801.(qfld).date;
    klgf1.(dqfld).data(end+1:end+n,1) = stn.KLB0801.(qfld).data;


    mrtf1.(srfld) = stn.k3mar101.(rfld);
    n = numel(stn.MKT0801.(rfld).date);
    mrtf1.(srfld).date(end+1:end+n,1) = stn.MKT0801.(rfld).date;
    mrtf1.(srfld).data(end+1:end+n,1) = stn.MKT0801.(rfld).data;
    % QA'd 10-minute data
    mrtf1.(sqfld) = stn.k3mar101.(qfld);
    n = numel(stn.MKT0801.(qfld).date);
    mrtf1.(sqfld).date(end+1:end+n,1) = stn.MKT0801.(qfld).date;
    mrtf1.(sqfld).data(end+1:end+n,1) = stn.MKT0801.(qfld).data;

    mrtf1.(drfld) = stn.k4mar101.(rfld);
    n = numel(stn.MKB0801.(rfld).date);
    mrtf1.(drfld).date(end+1:end+n,1) = stn.MKB0801.(rfld).date;
    mrtf1.(drfld).data(end+1:end+n,1) = stn.MKB0801.(rfld).data;
    % QA'd 10-minute data
    mrtf1.(dqfld) = stn.k4mar101.(qfld);
    n = numel(stn.MKB0801.(qfld).date);
    mrtf1.(dqfld).date(end+1:end+n,1) = stn.MKB0801.(qfld).date;
    mrtf1.(dqfld).data(end+1:end+n,1) = stn.MKB0801.(qfld).data;


    % QA'd interpolated hourly time-series
    klgf1.(sfld) = interp_ts(klgf1.(sqfld));
    klgf1.(dfld) = interp_ts(klgf1.(dqfld));
    mrtf1.(sfld) = interp_ts(mrtf1.(sqfld));
    mrtf1.(dfld) = interp_ts(mrtf1.(dqfld));

  end; %for cfld={'u','v','w','seapres','seatemp','xshore','lshore'};

  stn=[]; clear stn;

  [klgf1.cm_shallow_speed_raw,klgf1.cm_shallow_dir_raw] = ...
      ts_uv_to_spddir(klgf1.cm_shallow_u_raw,klgf1.cm_shallow_v_raw);
  [klgf1.cm_deep_speed_raw,klgf1.cm_deep_dir_raw] = ...
      ts_uv_to_spddir(klgf1.cm_deep_u_raw,klgf1.cm_deep_v_raw);

  [klgf1.cm_shallow_speed_10min,klgf1.cm_shallow_dir_10min] = ...
      ts_uv_to_spddir(klgf1.cm_shallow_u_10min,klgf1.cm_shallow_v_10min);
  [klgf1.cm_deep_speed_10min,klgf1.cm_deep_dir_10min] = ...
      ts_uv_to_spddir(klgf1.cm_deep_u_10min,klgf1.cm_deep_v_10min);

  [klgf1.cm_shallow_speed,klgf1.cm_shallow_dir] = ...
      ts_uv_to_spddir(klgf1.cm_shallow_u,klgf1.cm_shallow_v);
  [klgf1.cm_deep_speed,klgf1.cm_deep_dir] = ...
      ts_uv_to_spddir(klgf1.cm_deep_u,klgf1.cm_deep_v);


  [mrtf1.cm_shallow_speed_raw,mrtf1.cm_shallow_dir_raw] = ...
      ts_uv_to_spddir(mrtf1.cm_shallow_u_raw,mrtf1.cm_shallow_v_raw);
  [mrtf1.cm_deep_speed_raw,mrtf1.cm_deep_dir_raw] = ...
      ts_uv_to_spddir(mrtf1.cm_deep_u_raw,mrtf1.cm_deep_v_raw);

  [mrtf1.cm_shallow_speed_10min,mrtf1.cm_shallow_dir_10min] = ...
      ts_uv_to_spddir(mrtf1.cm_shallow_u_10min,mrtf1.cm_shallow_v_10min);
  [mrtf1.cm_deep_speed_10min,mrtf1.cm_deep_dir_10min] = ...
      ts_uv_to_spddir(mrtf1.cm_deep_u_10min,mrtf1.cm_deep_v_10min);

  [mrtf1.cm_shallow_speed,mrtf1.cm_shallow_dir] = ...
      ts_uv_to_spddir(mrtf1.cm_shallow_u,mrtf1.cm_shallow_v);
  [mrtf1.cm_deep_speed,mrtf1.cm_deep_dir] = ...
      ts_uv_to_spddir(mrtf1.cm_deep_u,mrtf1.cm_deep_v);

return;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INTERNAL FUNCTIONS

function stn = an_ncore_qa(stn)
% Apply simple Quality Assurance to Current Meter mooring sea pressure data
% in struct STN, returning STN with new "_qc" fields for u,v,w,seatemp,seapres.

  baddt = unique(floor(stn.cm_seapres_raw.date(abs(stn.cm_seapres_raw.data-median(stn.cm_seapres_raw.data))>5)));
  stn.cm_u_qc = subset_ts(stn.cm_u_raw,@(x)(find(~ismember(floor(x.date),baddt))));
  stn.cm_v_qc = subset_ts(stn.cm_v_raw,@(x)(find(~ismember(floor(x.date),baddt))));
  stn.cm_w_qc = subset_ts(stn.cm_w_raw,@(x)(find(~ismember(floor(x.date),baddt))));
  stn.cm_seatemp_qc = subset_ts(stn.cm_seatemp_raw,@(x)(find(~ismember(floor(x.date),baddt))));
  stn.cm_seapres_qc = subset_ts(stn.cm_seapres_raw,@(x)(find(~ismember(floor(x.date),baddt))));

return;
