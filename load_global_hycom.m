function station = load_global_hycom(station_or_stnm)
%function station = load_global_hycom(station_or_stnm)
%
% Load Global HYCOM analysis data from pre-saved MAT file into struct STATION.
%
% Last Saved Time-stamp: <Sun 2010-07-25 15:47:24 Eastern Daylight Time gramer>

  datapath = get_thesis_path('../data');

  if ( ischar(station_or_stnm) )
    stnm = station_or_stnm;
    station.station_name = stnm;
  elseif ( isstruct(station_or_stnm) )
    station = station_or_stnm;
    stnm = station.station_name;
  else
    error('First arg must be a station struct or station name string!');
  end;

  hycomfname = fullfile(datapath,[stnm '_hycom.mat']);
  rawstokesu = 'stokes_drift_u';
  if ( exist(hycomfname,'file') )
    disp(['Loading from presaved ' hycomfname]);
    load(hycomfname);

    dts = stn.global_hycom_u.date(1):(1/24):stn.global_hycom_u.date(end);
    station.global_hycom_u.date = dts';
    station.global_hycom_u.data = interp1(stn.global_hycom_u.date,...
                                          stn.global_hycom_u.data,...
                                          dts,'pchip',nan)';
    station.global_hycom_v.date = dts';
    station.global_hycom_v.data = interp1(stn.global_hycom_v.date,...
                                          stn.global_hycom_v.data,...
                                          dts,'pchip',nan)';
    stn = []; clear stn;
  elseif ( isfield(station,rawstokesu) )
    warning('Substituting all-zero surface currents for missing "%s"!', hycomfname);
    station.global_hycom_u.date = station.(rawstokesu).date;
    station.global_hycom_u.data = repmat(0,size(station.global_hycom_u.date));
    station.global_hycom_v = station.global_hycom_u;
  else
    warning('Global HYCOM current fields empty: missing "%s"!', hycomfname);
    station.global_hycom_u = struct('date',[],'data',[]);
    station.global_hycom_v = struct('date',[],'data',[]);
  end;

return;
