function stn = station_calc_sgs_diffusion(stn,sfld,dtfld,tffld,lplfld,sgsfld,kfld,bdgtfld)
%function stn = station_calc_sgs_diffusion(stn,sfld,dtfld,tffld,lplfld,sgsfld,kfld,bdgtfld)
%
% Calculate "sub-grid-scale" heat diffusion required to balance heat budget.
% BDGTFLD arg is optional; if present, create a final, balanced budget field. 
%
% SAMPLE CALL:
% >> stn = station_calc_sgs_diffusion(stn,'ndbc_sea_t',...
%                                     'fkeys_hycom_qelnetqf', ...
%                                     'fkeys_hycom_seatemp_field', ...
%                                     'fkeys_hycom_seatemp_laplacian', ...
%                                     'fkeys_hycom_sgs_diffused_heat', ...
%                                     'fkeys_hycom_sgs_thermal_diffusivity', ...
%                                     'fkeys_hycom_nn30a_final_budget');
%
% Last Saved Time-stamp: <Wed 2012-03-21 12:49:40  Lew.Gramer>

  if ( ~exist('bdgtfld') || isempty(bdgtfld) )
    bdgtfld = [];
  end;

  if ( ~isfield(stn,lplfld) )
    warning('Midpoint-interpolating Laplacian');
    % Spline-fit the central pixel of the six-hourly FKEYS (or daily GoM)
    % sea temperature field Laplacian into an hourly time series 
    midx = round(size(stn.(tffld).laplacian,2)/2);
    midy = round(size(stn.(tffld).laplacian,3)/2);
    rawdts = stn.(tffld).date;
    rawdat = squeeze(stn.(tffld).laplacian(:,midx,midy));

    dts = [rawdts(1):(1/24):rawdts(end)]';
    dat = spline(rawdts,rawdat,dts);
    stn.(lplfld).date = dts;
    stn.(lplfld).data = dat;
  end;

  % Calculate an "effective" thermal diffusivity based on the residual from
  % our total heat budget, and actual hourly sea temperature change

  % DH = T(i+1) - T(i) - dTdt	Sub-gridscale (effective) Diffused Heat [K/hr]
  % DH/3600 = KHs * del2(T)	[K/s] = [m^2/s] * [K/m^2]
  % KHs = DH/(3600*del2(T))	[m^2/s] = [K/hr]/[s*K/hr*m^2]

  [ix1,ix2] = intersect_dates(stn.(sfld).date,stn.(dtfld).date);
  p = stn.(sfld).data(ix1)+stn.(dtfld).data(ix2);
  stn.(sgsfld).date = stn.(sfld).date(ix1(1:end-1));
  stn.(sgsfld).data = stn.(sfld).data(ix1(2:end))-p(1:end-1);
  % Handle time series gaps in sea temperature and heat budget
  badix = find( (diff(stn.(sfld).date(ix1)) > (2/24)) | ...
                (diff(stn.(dtfld).date(ix2)) > (2/24)) );
  stn.(sgsfld).date(badix) = [];
  stn.(sgsfld).data(badix) = [];

  [ix1,ix2] = intersect_dates(stn.(sgsfld).date,stn.(lplfld).date);
  stn.(kfld).date = stn.(sgsfld).date(ix1);
  stn.(kfld).data = stn.(sgsfld).data(ix1)./(stn.(lplfld).data(ix2).*3600);

  % Handle time series gaps in SGS diffused heat, and spline-fit hourly Laplacian
  badix = find( (diff(stn.(sgsfld).date(ix1)) > (2/24)) | ...
                (diff(stn.(lplfld).date(ix2)) > (2/24)) );
  stn.(kfld).date(badix) = [];
  stn.(kfld).data(badix) = [];


  % Remove unphysical estimates of K_heat_SGS
  badix = find(stn.(kfld).data < 0);
  %DEBUG:  disp(['NOT Removing ' num2str(numel(badix)) ' unphysical K_H estimates']);
  disp(['Removing ' num2str(numel(badix)) ' unphysical K_H estimates']);
  stn.(kfld).date(badix) = [];
  stn.(kfld).data(badix) = [];


  % Remove excessively energetic estimates
  % Assume KHs ~ U_SGS*L_SGS: U_SGS <= KHs/1km for FKEYS and USF AVHRR, 4km for GoM
  dx = sw_dist(stn.(tffld).lat([1 1]),stn.(tffld).lon([1 2]),'km').*1e3;
  Usgs = stn.(kfld).data./dx;
  badix = find(Usgs > 1.0);
  %DEBUG:
  disp(['NOT Removing ' num2str(numel(badix)) ' excessive K_H estimates']);
  % disp(['Removing ' num2str(numel(badix)) ' excessive K_H estimates']);
  % stn.(kfld).date(badix) = [];
  % stn.(kfld).data(badix) = [];

  if ( ~isempty(bdgtfld) )
    [ix1,ix2] = intersect_dates(stn.(dtfld).date,stn.(sgsfld).date);
    stn.(bdgtfld).date = stn.(dtfld).date(ix1);
    stn.(bdgtfld).data = stn.(dtfld).data(ix1) + stn.(sgsfld).data(ix2);
  end;

return;
