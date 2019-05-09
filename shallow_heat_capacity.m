function stn = shallow_heat_capacity(stn,topfld,btmfld,zfld,q0fld)

%   'bic_surf_par'
%   'bic_shallow_par'
%   'tmd_tide_i_depth'

  stn = calc_kd(stn,topfld,btmfld,zfld);

  % PAR attenuation coefficient in water column
  Kd=0.06;
  kdfld = ['kd_' topfld '_' btmfld];
  if ( isfield(stn,kdfld) )
    Kd = nanmean(real(stn.(kdfld).data));
    Kd(Kd <= 0) = [];
  else
    Kd = 0.04;
  end;

  % Solar zenith angle
  theta=10;

  % Sea-bottom mean reflectivity
  Ab=0.4;

  q0 = 400;

  sun_angle_correction = 1.0;
  if ( isfield(stn,'lon') )
    [yds,yrs] = get_yearday(stn.(newfld).date);
    [theta,ig] = soradna1(yds,yrs,-stn.lon,stn.lat);
    sun_angle_correction = secd(90-theta(ix));
  end;

  % Assume 50% of insolation in NIR, total absorbed in the water column
  QSWI = (q0/2) + ...
         (q0/2)*(1 - exp(-Kd*z*sun_angle_correction) + ...
                 (Ab*exp(-Kd*z*sun_angle_correction)));

return;
