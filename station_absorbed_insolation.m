function stn = station_absorbed_insolation(stn,aifld,qswfld,zfld,topfld,btmfld,gamfld,qbfld,opts,diagfld)
%function stn = station_absorbed_insolation(stn,aifld,qswfld,zfld,topfld,btmfld,gamfld,qbfld,opts)
%
% Calculate insolation absorbed by the ocean, STN.(AIFLD) [W/m2], from net sea
% surface insolation STN.(QSWFLD) [W/m2], depth STN.(ZFLD) [m], and either a
% constant or seasonal Kd_PAR; or (if TOPFLD and BTMFLD both given and valid)
% a Kd_PAR calculated from the upper and lower PAR data in STN.(TOPFLD) and
% STN.(BTMFLD) fields, resp., and depth field STN.(zfld).  Assumes bottom
% reflectance for PAR based on constant reflectances for sand/rubble and all
% other bottom types, resp.  Tries calling STATION_SAND_COVER to get fraction
% (0-1) of site sea-floor that is sand/rubble; otherwise assumes sand=60%. If
% optional GAMFLD string given, calculated time series of absorption factor
% (gamma) is returned in new field STN.(GAMFLD). If QBFLD given, insolation
% absorbed by the benthos is returned in new time series field STN.(QBFLD).
%
% Optionally, ZFLD may also be a scalar, specifying a constant site depth.
%
% Optional OPTS struct may specify any or all of the following options:
%  OPTS.par_per_insol : fraction (0-1) of incident short-wave radiation taken
%        to be in PAR band (400-700nm). DEFAULT is 0.501.
%  OPTS.uv_per_insol : fraction (0-1) of incident short-wave radiation taken
%        to be in near-UV (NUV=UVA+UVB) band (280-400nm). DEFAULT is 0.048.
%  OPTS.Ppen : fraction (0-1) of incident short-wave rad. in PAR+NUV (i.e.,
%        penetrative). DEFAULT: OPTS.par_per_insol+OPTS.uv_per_insol=0.549.
%  OPTS.kd : If scalar, a constant Kd (e-folding depth) for PAR and near-UV.
%        If vector length==2 (or 3 or 4), min and max range (and peak day,
%        period in days) for a sine wave-like time series of Kd. DEFAULT
%        is annual variability with peak year-day 91.5 (12 UT April 01). See
%        BUILD_CLIM_OPT for all the option formats acceptable for OPTS.kd. 
%  OPTS.allow_chlor_kd : if true *and* field STN.amodis_chlor_a is present,
%        then attempt to use methodology of Ohlmann (2003) to calculate
%        Kd time series from satellite chlorophyll /a/ concentrations.
%  OPTS.sand_fraction : fraction (0-1) of site which should be considered
%        sand (i.e., highly reflective) as opposed to all other benthic
%        habitat types (mud, hard-bottom, live coral, seagrass).
%  OPTS.bottom_reflectance : a reflectance Ab (0-1) for sea floor overriding
%        that calculated according to sand fraction (see above).
%  OPTS.absorption_factor_gamma : an absorption factor for net sea-surface
%        short-wave radiation which overrides that calculated based on
%        attenuation coefficient Kd or sea-floor reflectance Ab (above).
%
% NOTE: This function does not account for benthic *turbulent* or *long-wave*
% heat exchange. Calling it on a shallow site without calculating these terms
% separately may cause "frozen reef syndrome" in your heat budget: BEWARE!
%
% Suggested inputs for various SEAKEYS and ICON stations:
%  ICON: zfld='ctd_deep_i_depth' and 'kd_bic_surf_par_bic_deep_par'
%  MLRF1: zfld='tmd_tide_i_depth' and 'kd_bic_surf_par_bic_shallow_par'
%  SMKF1: zfld='tide_i_depth' and Kd=constant 0.30 or seasonal 0.15-0.45
%  LONF1: zfld='tide_i_depth' and Kd=constant 0.30 or seasonal 0.15-0.45
%  Everywhere else: zfld='tmd_tide_i_depth' and Kd=constant or seasonal
%
% CALLS: NANMEAN (Statistics Toolbox); GET_OPT; CALC_KD, STATION_SAND_COVER,
%   INTERSECT_DATES (Ecoforecasts Toolbox); SORADNA1 (Air_Sea Toolbox).
%
% Last Saved Time-stamp: <Tue 2012-11-27 12:44:20 Eastern Standard Time lew.gramer>


% Additional private arg DIAGFLD names another new sub-struct returned in STN
% with diagnostics: sun altitude and correction, tau (one-way attenuation).

  if ( ~exist('opts','var') || isempty(opts) )
    opts = [];
  end;

  doDebug = get_opt(opts,'kd_debug',true);

  % Conversion fraction from Insolation (global radiation) to PAR [W/m^2]
  % % Original default 0.473 was from Papaioannou, Papanikolaou and Retalis (1993)
  % PAR_PER_INSOL = get_opt(opts,'par_per_insol',0.473);
  % Hourly mean for overcast skies, per Jacovides et al., 2003
  PAR_PER_INSOL = get_opt(opts,'par_per_insol',0.501);

  % % Per Canada et al., 2003
  % UV_PER_INSOL = get_opt(opts,'uv_per_insol',0.050);
  % % Mean per Ecobedo et al., 2009 +/- 18%
  % UV_PER_INSOL = get_opt(opts,'uv_per_insol',0.049);
  % Mean fraction daily dose UV(A+B) of daily dose insolation, for the annual
  % range of Keys total cloud cover given by ERAI, per Leal et al., 2011
  UV_PER_INSOL = get_opt(opts,'uv_per_insol',0.048);

  % Ppen = Ppar + Pnuv
  Ppen = get_opt(opts,'Ppen',PAR_PER_INSOL+UV_PER_INSOL);
  if ( isempty(Ppen) )
    Ppen = PAR_PER_INSOL+UV_PER_INSOL;
  end;

  if ( ~exist('topfld','var') || isempty(topfld) || ...
       ~exist('btmfld','var') || isempty(btmfld) )
    kdfld = 'kd_field_does_not_exist';
  else
    stn = calc_kd(stn,topfld,btmfld,zfld);
    kdfld = ['kd_' topfld '_' btmfld];
  end;

  % What does caller want us to use for station depth?
  if ( ischar(zfld) )
    [zix,qswix] = intersect_dates(stn.(zfld).date,stn.(qswfld).date);
    z = stn.(zfld).data(zix);
  else
    qswix = 1:length(stn.(qswfld).date);
    if ( isfield(stn,'depth') )
      z = stn.depth;
    else
      z = 1.0;
    end;
  end;

  % PAR attenuation coefficient in water column
  Kd = get_opt(opts,'kd',[]);
  allowChlorKd = get_opt(opts,'allow_chlor_kd',false);

  if ( isempty(Kd) )
    if ( isfield(stn,kdfld) )
      % Calculated from in situ light data
      Kds = real(stn.(kdfld).data);
      Kds(1 > Kds | Kds <= 0) = [];
      Kd = nanmean(Kds);

    elseif ( allowChlorKd && isfield(stn,'amodis_chlor_a') )
      warning('This algorithm is in desparate need of debugging!');
      %% Ohlmann [2003] parameterization of Kd from chlorophyll a [mg/m^3]
      %% See Sweeney et al. [2005] for field comparison results
      stn = verify_variable(stn,'amodis_chlor_a_1_hour_spline');
      chlz = z;
      chl = stn.amodis_chlor_a_1_hour_spline.data./10;
      chldts = stn.amodis_chlor_a_1_hour_spline.date;
      if ( ischar(zfld) )
        [chlzix,chlix] = intersect_dates(stn.(zfld).date(zix),stn.amodis_chlor_a_1_hour_spline.date);
        chlz = z(chlzix);
        chl = chl(chlix);
        chldts = chldts(chlix);
      end;
      zeta1=0.015+(0.176*((0.462*chl).^0.5));
      zeta2=0.688+(0.060.*log(0.125*chl));
      V1=0.571 + (0.025.*log(0.149*chl));
      V2=0.22374007 - (0.01609314.*log(2.33384554*chl));
      Ifac = (V1.*exp(-chlz./zeta1)) + (V2.*exp(-chlz./zeta2));
      stn.kd_amodis_chlor_a.date=chldts;
      stn.kd_amodis_chlor_a.data=-log(Ifac)./chlz;
      [kdix,ig] = intersect_dates(stn.kd_amodis_chlor_a.date,stn.(qswfld).date);
      Kd = stn.kd_amodis_chlor_a.data(kdix);
      if ( length(Kd) < length(qswix) )
        warning('Insufficient time series of Aqua MODIS chlor_a: using median.');
        Kd = nanmean(Kd);
      end;

    else
      % Use default Kd that minimizes error for *all* SEAKEYS sea temperatures 
      if ( doDebug )
        disp('Using constant Kd = 0.30');
      end;
      Kd=0.30;
    end;

  else %if ( isempty(Kd) ) else

    if ( ischar(Kd) )
      % Mean calculated from in situ light data
      Kds = real(stn.(Kd).data);
      Kds(1 > Kds | Kds <= 0) = [];
      Kd = nanmean(Kds);

    elseif ( isnumeric(Kd) )
      if ( numel(Kd) == 2 )
        % DEFAULT peak year-day for seasonal climatology is in MARCH
        Kd(3) = 91;
      end;
      Kd = build_clim_opt(Kd,'Kd',stn.(qswfld).date(qswix),doDebug);

    elseif ( is_ts(Kd) || ( iscell(Kd) && numel(Kd)==2 ))
      Kd = build_clim_opt(Kd,'Kd',stn.(qswfld).date(qswix),doDebug);

    else
      error('OPTS.Kd may be a time series field, name, cell ts-function-pair, or numeric scalar/vector!');

    end; %if ( ischar(Kd) ) elseif ( isnumeric(Kd) ) else

  end; %if ( isempty(Kd) ) else


  % Sea-floor habitat types
  frac = get_opt(opts,'sand_fraction',[]);
  if ( isempty(frac) )
    try
      frac = station_sand_cover(stn);
    catch
      % Default 60% sand cover, based on Keys-wide estimates of hard-bottom/sea
      % grass mean cover >=30% and living reef almost universally down to <10%
      frac = 0.60;
    end;
  end;

  % Sea-bottom mean reflectance in VISIBLE/UV light:
  % FWYF1 0.235, MLRF1 0.235, LONF1 0.169, SMKF1 0.235, LOOE1 0.400
  Ab = get_opt(opts,'bottom_reflectance',[]);
  if ( isempty(Ab) )
    % DEFAULT reflectances for sand and non-sand from E. Hochberg (pers
    % comm). See also Louchard et al (2003): Oolitic-skeletal sediment ~0.50;
    % T. testudinum ~0.08. Zhang Voss Reid Louchard (2003) find bidirectional
    % coeff. "REFF" for ooid sand ~0.46, algal film over sand ~0.25.
    Ab = (0.40 * frac) + (0.07 * (1-frac));
  end;

  dts = stn.(qswfld).date(qswix);
  qsw = stn.(qswfld).data(qswix);

  % Solar zenith angle correction. (Note - this is necessary and distinct
  % from the near-identical calculation which is done in CALC_KD above.)
  sun_alt = repmat(90.0,size(dts));
  if ( isfield(stn,'lon') )
    % Center hourly means on the half-hour
    [yds,yrs] = get_yearday(dts-(0.5/24));
    [sun_alt,ig] = soradna1(yds,yrs,-stn.lon,stn.lat);
  end;
  sun_alt_correction = secd(90-sun_alt);

  % Fraction of incident light (PAR) *NOT* absorbed by water column in each direction
  tau = exp(-Kd .* z .* sun_alt_correction);
  % For light scattered by benthos - assumes a "mean beam angle" of 45o
  tau_scatter = exp(-Kd .* z .* secd(90-45));

  % Basic sanity check
  tau(0>tau | tau>1) = 0;
  tau_scatter(0>tau_scatter | tau_scatter>1) = 0;

  %DEBUG:  fmg; plot(dts,tau); title('\tau');

  % The following differs significantly from, yet benefitted from a reading
  % of Warrior and Carder (2007) and citations there, e.g., Lee et al (2001)

% ??? VERIFY against Lewis, Wei, van Dommelen, Voss (2011)

  % % Assume <50% of insolation in near infrared (NIR), i.e., totally absorbed.
  % % Calculate (1 - Kd_PAR) absorption of the other <50%, in both directions.
  % % Assumes SPECULAR REFLECTION of insolation by the sea-floor; also assumes
  % % any insolation not reflected by the sea-floor is immediately absorbed and
  % % conducted into the sea-bed FOREVER TO DISAPPEAR: therefore, a separate
  % % calculation of benthic-water column heat exchange *must* be done.
  % gam = ( 1 - Ppen ) + ( Ppen.*(1 - tau + (tau.*Ab.*(1-tau_scatter))) );

  gam = get_opt(opts,'absorption_factor_gamma',[]);
  if ( isempty(gam) )
    % NIR component (fully absorbed)
    gam = (1 - Ppen);
    % Downward PAR and NUV absorbed
    gam = gam + (Ppen.*(1 - tau));

%%%%DEBUG???
    % % PAR reflected and scattered from benthos (assume ALL absorbed)
    % gam = gam + (Ppen.*tau.*Ab);
    % % PAR reflected specularly from benthos and absorbed (reuse TAU)
    % gam = gam + (Ppen.*tau.*Ab.*(1-tau));
    % PAR scattered from benthos and absorbed
    gam = gam + (Ppen.*tau.*Ab.*(1-tau_scatter));
  end;

  % Basic sanity check
  gam(eps>=gam | gam>1) = 1;

  stn.(aifld).date = dts;
  stn.(aifld).data = gam .* qsw;

  % When the sun is too low for reasonable estimates, ignore insolation
  badix = find(sun_alt <= -6);	% Normal rise/set rate 15 degrees/hour
  stn.(aifld).data(badix) = 0;

  % If caller is interested in absorption factor for future reference
  if ( exist('gamfld','var') && ~isempty(gamfld) )
    if ( isscalar(gam) )
      gam = repmat(gam,size(dts));
    end;
    stn.(gamfld).date = dts;
    stn.(gamfld).data = gam;
    stn.(gamfld).data(badix) = 1;
  end;

  % If caller is interested in insolation absorbed by sea bed. (Assumes all
  % NIR radiation is absorbed in water column: only PAR/NUV reaches bottom.)
  if ( exist('qbfld','var') && ~isempty(qbfld) )
    stn.(qbfld).date = dts;
    stn.(qbfld).data = qsw .* Ppen .* tau .* (1 - Ab);
    stn.(qbfld).data(badix) = 0;
  end;

  % DEBUG: If caller asks for diagnostic time series, return them
  if ( exist('diagfld','var') && ~isempty(diagfld) )
    stn.(diagfld).Ppen = Ppen;
    stn.(diagfld).sun_alt.date = dts;
    stn.(diagfld).sun_alt.data = sun_alt;
    stn.(diagfld).sun_alt_correction.date = dts;
    stn.(diagfld).sun_alt_correction.data = sun_alt_correction;
    stn.(diagfld).sun_alt_correction.data(badix) = 0;
    stn.(diagfld).tau.date = dts;
    stn.(diagfld).tau.data = tau;
    stn.(diagfld).tau.data(badix) = 0;
    stn.(diagfld).tau_scatter.date = dts;
    stn.(diagfld).tau_scatter.data = tau_scatter;
    stn.(diagfld).tau_scatter.data(badix) = 0;
    stn.(diagfld).Kd = stn.(diagfld).tau;
    stn.(diagfld).Kd.data(:) = Kd;
    % Shortwave radiation reaching the benthos (whether reflected or absorbed)
    stn.(diagfld).bottom_dsr.date = dts;
    stn.(diagfld).bottom_dsr.data = qsw .* Ppen .* tau;
    stn.(diagfld).bottom_dsr.data(badix) = 0;
    % Shortwave radiation lost back into the air
    stn.(diagfld).lost_sr.date = dts;
    stn.(diagfld).lost_sr.data = (1-gam) .* qsw;
  end;

return;
