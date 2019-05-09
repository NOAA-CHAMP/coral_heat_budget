function stn = fix_gom_hycom_field(stnm)
%function stn = fix_gom_hycom_field(stnm)
%
% Add lat and lon fields for Gulf of Mexico 1/25-degree HYCOM sea temperature
% field, rename '.data' field as '.field' field. Calculate 'speed' and 'dir'
% time series fields from STN.gom_hycom_u and STN.gom_hycom_v time series.

error('This is a "run-once" fix-it function! Do not run again...');

  datapath = get_thesis_path('../data');

  stn = get_gom_hycom(stnm);

  % Calculate speed and direction from model U and V currents
  if ( isfield(stn,'gom_hycom_u') && isfield(stn,'gom_hycom_v') )
    stn.gom_hycom_speed.date = stn.gom_hycom_u.date;
    stn.gom_hycom_speed.data = uv_to_spd(stn.gom_hycom_u.data,stn.gom_hycom_v.data);
    stn.gom_hycom_dir.date = stn.gom_hycom_u.date;
    stn.gom_hycom_dir.data = uv_to_dir_curr(stn.gom_hycom_u.data,stn.gom_hycom_v.data);
  end;

  if ( ~exist('fldnm','var') || isempty(fldnm) )
    fldnm = 'gom_hycom_seatemp_field';
  end;

  % Modify STN.gom_hycom_seatemp_field struct: change field .data to
  % .field, and add Nx1 fields .lon and .lat, before saving to MAT.
  stn = query_gom_hycom_station_field_coords(stn);

  matfname = fullfile(datapath,[lower(stn.station_name) '_gom_hycom.mat']);
  disp(['RE-saving to MAT file ' matfname]);
  save(matfname,'stn');

return;
