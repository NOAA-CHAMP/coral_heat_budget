function fix_erai_station(stnm)

  error('%s: This fix-it script only needs to be run ONCE per broken station!',stnm);

  datapath = get_thesis_path('../data');

  matfname = fullfile(datapath,[stnm '_erai.mat']);
  if ( ~exist(matfname,'file') )
    error('MAT file not found! "%s"',matfname);
  end;

  load(matfname,'station');

  station.raw_erai_uvW = station.raw_erai_durf;
  station.erai_uvW = station.erai_durf;
  station = rmfield(station,{'erai_durf','raw_erai_durf'});

  % Calculate some additional flux fields from those given
  station.raw_erai_usrf.date = station.raw_erai_dsrf.date;
  station.raw_erai_usrf.data = ...
      station.raw_erai_dsrf.data - station.raw_erai_srf.data;
  station.raw_erai_ulrf.date = station.raw_erai_dlrf.date;
  station.raw_erai_ulrf.data = ...
      station.raw_erai_dlrf.data - station.raw_erai_lrf.data;

  rawflds = {'raw_erai_usrf','raw_erai_ulrf',};
  for fldix = 1:length(rawflds)
    rawfld = rawflds{fldix};
    fld = rawfld(5:end);

    rawdts = station.(rawfld).date;
    dts = [rawdts(1):(1/24):rawdts(end)]';
    ndts = length(dts);

    rawdat = station.(rawfld).data;
    station.(fld).date(1:ndts,1) = dts(:);
    station.(fld).data(1:ndts,1) = interp1(rawdts,rawdat,dts(:),'spline');
    station = filter_gaps(station,rawfld,fld,[],(6/24),[],nan);

    station.(fld).data(station.(fld).data<0) = 0;
  end;


  % % Very simplistic QA:
  % % Certain spline-interpolated fields should never be negative
  % station.erai_cloud_cover.data(station.erai_cloud_cover.data<0) = 0;
  % station.erai_conv_precip.data(station.erai_conv_precip.data<0) = 0;
  % station.erai_durf.data(station.erai_durf.data<0) = 0;
  % % (Evaporation should never be positive!)
  % station.erai_evap.data(station.erai_evap.data>0) = 0;
  % station.erai_parW.data(station.erai_parW.data<0) = 0;
  % station.erai_runoff.data(station.erai_runoff.data<0) = 0;
  % station.erai_sunshine_duration.data(station.erai_sunshine_duration.data<0) = 0;
  % station.erai_clear_sky_srf.data(station.erai_clear_sky_srf.data<0) = 0;
  % station.erai_dsrf.data(station.erai_dsrf.data<0) = 0;
  % station.erai_srf.data(station.erai_srf.data<0) = 0;
  % station.erai_dlrf.data(station.erai_dlrf.data<0) = 0;
  % station.erai_precip.data(station.erai_precip.data<0) = 0;

  % % Additional calculations for some interpolated fields
  % % Dew P -> Rel Humid -> Spec Humid
  % station.erai_relhumid.date = station.erai_dew_t.date;
  % station.erai_relhumid.data = ...
  %     dewp_to_relhumid(station.erai_air_t.data,station.erai_dew_t.data);
  % station.erai_spechumid.date = station.erai_relhumid.date;
  % station.erai_spechumid.data = ...
  %     relhumid_to_spechumid(station.erai_air_t.data,station.erai_relhumid.data);

  disp('Resaving MAT file...');
  save(matfname,'station');
  station = []; clear station;

return;
