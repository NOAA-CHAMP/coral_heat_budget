function [fkeys_hycom_lons,fkeys_hycom_lats] = save_fkeys_hycom_coords

  warning('Should only be run ONE TIME per installation!');
  pause;

  datapath = get_thesis_path('../data');
  avhrrpath = fullfile(datapath,'hycom','FKEYS');

  fname = fullfile(avhrrpath,'303_archv.2004_001_06_3zu.nc');
  nc = mDataset(fname);
  fkeys_hycom_lons = cast(nc{'Longitude'}(:,:),'double');
  fkeys_hycom_lats = cast(nc{'Latitude'}(:,:),'double');
  close(nc); clear nc;

  save(fullfile(datapath,'fkeys_hycom_coords.mat'),'fkeys_hycom_lons','fkeys_hycom_lats');

return;
