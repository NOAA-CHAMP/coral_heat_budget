function stn = calc_sgs(stn)
%function stn = calc_sgs(stn)
%
% Calculate "sub-grid-scale" heat diffusion required to balance heat budget.
%
% Last Saved Time-stamp: <Wed 2010-12-15 16:11:00 Eastern Standard Time gramer>

  % Spline-fit six-hourly FKEYS sea temperature Laplacian to hourly time series 
  rawdts = stn.fkeys_hycom_seatemp_field.date;
  rawdat = squeeze(stn.fkeys_hycom_seatemp_field.laplacian(:,9,9));

  dts = [rawdts(1):(1/24):rawdts(end)]';
  dat = spline(rawdts,rawdat,dts);
  stn.fkeys_hycom_seatemp_laplacian.date = dts;
  stn.fkeys_hycom_seatemp_laplacian.data = dat;

  % Calculate an "effective" thermal diffusivity based on the residual from
  % our total heat budget, and actual hourly sea temperature change

  % DH = T(i+1) - T(i) - dTdt	Sub-gridscale (effective) Diffused Heat [K/hr]
  % DH/3600 = KHs * del2(T)	[K/s] = [m^2/s] * [K/m^2]
  % KHs = DH/(3600*del2(T))	[m^2/s] = [K/hr]/[s*K/hr*m^2]

  [ix1,ix2] = intersect_dates(stn.ndbc_sea_t.date,stn.fkeys_hycom_qelnetqf.date);
  p = stn.ndbc_sea_t.data(ix1)+stn.fkeys_hycom_qelnetqf.data(ix2);
  stn.fkeys_hycom_sgs_diffused_heat.date = stn.ndbc_sea_t.date(ix1(1:end-1));
  stn.fkeys_hycom_sgs_diffused_heat.data = stn.ndbc_sea_t.data(ix1(2:end))-p(1:end-1);
  % Handle time series gaps in sea temperature and heat budget
  badix = find( (diff(stn.ndbc_sea_t.date(ix1)) > (2/24)) | ...
                (diff(stn.fkeys_hycom_qelnetqf.date(ix2)) > (2/24)) );
  stn.fkeys_hycom_sgs_diffused_heat.date(badix) = [];
  stn.fkeys_hycom_sgs_diffused_heat.data(badix) = [];

  [ix1,ix2] = intersect_dates(stn.fkeys_hycom_sgs_diffused_heat.date,stn.fkeys_hycom_seatemp_laplacian.date);
  stn.fkeys_hycom_sgs_thermal_diffusivity.date = stn.fkeys_hycom_sgs_diffused_heat.date(ix1);
  stn.fkeys_hycom_sgs_thermal_diffusivity.data = stn.fkeys_hycom_sgs_diffused_heat.data(ix1)./(stn.fkeys_hycom_seatemp_laplacian.data(ix2).*3600);

  % Handle time series gaps in SGS diffused heat and spline-fit hourly Laplacian
  badix = find( (diff(stn.fkeys_hycom_sgs_diffused_heat.date(ix1)) > (2/24)) | ...
                (diff(stn.fkeys_hycom_seatemp_laplacian.date(ix2)) > (2/24)) );
  stn.fkeys_hycom_sgs_thermal_diffusivity.date(badix) = [];
  stn.fkeys_hycom_sgs_thermal_diffusivity.data(badix) = [];

  % Remove unphysical estimates
  badix = find(stn.fkeys_hycom_sgs_thermal_diffusivity.data <= eps);
  %DEBUG:
  disp(['Removing ' num2str(numel(badix)) ' unphysical K_H estimates']);
  stn.fkeys_hycom_sgs_thermal_diffusivity.date(badix) = [];
  stn.fkeys_hycom_sgs_thermal_diffusivity.data(badix) = [];


  % Remove excessively energetic estimates
  % Assume KHs ~ U_SGS*L_SGS, or U_SGS = KHs/1km
  Usgs = stn.fkeys_hycom_sgs_thermal_diffusivity.data/1e3;
  badix = find(Usgs > 1.0);
  %DEBUG:
  disp(['NOT Removing ' num2str(numel(badix)) ' excessive K_H estimates']);
%   stn.fkeys_hycom_sgs_thermal_diffusivity.date(badix) = [];
%   stn.fkeys_hycom_sgs_thermal_diffusivity.data(badix) = [];

  multiplot_station(stn,{'ndbc_sea_t_diff','fkeys_hycom_qelnetqf','fkeys_hycom_seatemp_laplacian','fkeys_hycom_sgs_diffused_heat','fkeys_hycom_sgs_thermal_diffusivity'});
  xlim([datenum(2004,12,9,6,0,0),datenum(2004,12,19,18,0,0)]);

return;
