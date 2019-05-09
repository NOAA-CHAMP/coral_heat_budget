function stn = repro_ms(stn)

  afld = 'ndbc_air_t';
  sfld = 'ndbc_sea_t';
  pfld = 'ndbc_barom';
  dfld = 'ndbc_dew_t';
  hfld = 'tmd_tide_i_depth';
  qafld = 'ndbc_spechumid';
  qsfld = 'ndbc_sea_spechumid';
  Wfld = 'ndbc_wind1_speed';
  ufld = 'fkeys_hycom_u';
  vfld = 'fkeys_hycom_v';
  Tfld = 'fkeys_hycom_seatemp_field';
  whfld = 'ndbc_air_t';
  wpfld = 'ndbc_air_t';
  wdfld = 'ndbc_air_t';

  dsrfld = 'erai_dsrf';
  usrfld = 'erai_usrf';
  dlrfld = 'erai_bulk_dlrf';
  ulrfld = 'erai_bulk_ulrf';

  qlhfld = 'ndbc_erai_30a_lhf';
  qshfld = 'ndbc_erai_30a_shf';
  qrhfld = 'ndbc_erai_30a_rhf';

  q0fld = 'ndbc_erai_30a_nhf';

  udTfld = 'fkeys_hycom_advected_heat';
  kd2Tfld = 'fkeys_hycom_diffused_heat';

return;
