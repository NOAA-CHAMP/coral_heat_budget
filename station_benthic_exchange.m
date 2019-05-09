function stn = station_benthic_exchange(stn,tfld,ufld,vfld,qbfld,tbfld,qbofld,opts,diagfld)
%function stn = station_benthic_exchange(stn,tfld,ufld,vfld,qbfld,tbfld,qbofld,opts)
%
% Calculate benthic temperature STN.(TBFLD) and heat flux into water column
% STN.(QBOFLD), from sea temperature STN.(TFLD), hourly currents UFLD,VFLD,
% and benthic absorbed short-wave radiation, STN.(QBFLD), assuming: initial
% equilibrium between benthic and water temperature; constant (mean) values
% for seabed density, heat capacity, long-wave emissivity, and conductivity.
%
% Optional OPTS struct may specify any of the following options:
%  OPTS.sand_fraction : fraction (0-1) of benthic active layer that is sand-
%      or coral rock-like (i.e., highly porous) vs. mud-like (low porosity).
%  OPTS.benthic_layer_depth : depth [m] of "thermally active" substrate
%  OPTS.benthic_boundary_layer_depth : depth [m] of substrate that we assume
%      to be in constant thermal equilibrium with water (DEFAULT 0.03).
%  OPTS.b_convective_coefficient : exchange coefficient for estimating the
%      convective heat flux from benthos to water column (DEFAULT 1.7e-4).
%
% Last Saved Time-stamp: <Sat 2012-03-31 02:26:37  Lew.Gramer>

% Additional private arg DIAGFLD names another new sub-struct returned in STN
% with diagnostics: all sub-components of benthos-water heat exchange model.

  %DEBUG:  tic,

  if ( ~exist('opts','var') || isempty(opts) )
    opts = [];
  end;

  doDebug = get_opt(opts,'benthic_debug',true);

  stn = verify_variable(stn,tfld);
  stn = verify_variable(stn,ufld);
  stn = verify_variable(stn,vfld);
  stn = verify_variable(stn,qbfld);

  [tix,uix,vix,qbix] = ...
      intersect_all_dates([],stn.(tfld).date,stn.(ufld).date,stn.(vfld).date,stn.(qbfld).date);

  dts = stn.(tfld).date(tix);
  t = stn.(tfld).data(tix);
  u = stn.(ufld).data(uix);
  v = stn.(vfld).data(vix);
  qb = stn.(qbfld).data(qbix);


  % Sea-bottom habitat types
  frac = get_opt(opts,'sand_fraction',[]);
  if ( isempty(frac) )
    try
      frac = station_sand_cover(stn);
    catch
      % Hard-bottom/sea grass mean >=30%, living reef down to <10%
      frac = 0.6;
    end;
  end;
  if ( doDebug )
    disp([mfilename ': Using sand fraction ' num2str(frac) '']);
  end;


  % Approximate seawater density [kg/m^3] and specific heat capacity [J/kg/K]
  rhow = 1.025e3;
  Cpw = 4e3;


  % NOTE: Many of the following assumptions were validated with literature
  % (as cited below), and with an inter-tidal mud study (Guarini et al 1997)

  % Benthic heating from incident insolation. ASSUMES benthic "active layer"
  % under sea bed is ~3m deep, and contains *a* mix of marine sediment and
  % coral rock: rock density 2600 [kg/m3; Hughes 1987] specific heat capacity
  % ~0.2 times that of seawater; marine sediment is 70% sandy clay, 30% water
  % mix by mass; dry sandy clay has density ~2600, specific heat ~0.3 times
  % seawater. Our rho_b, Cp_b are volume averages of the above properties.
  hb = get_opt(opts,'benthic_layer_depth',[]);
  if ( isempty(hb) )
    hb = 2;
  end;
  if ( doDebug )
    disp([mfilename ': Using thermal active depth ' num2str(hb) 'm']);
  end;

  % Areal-average substrate density [kg/m^3] and specific heat capacity [J/kg/K]
  rho_rock = 2.6.*rhow;
  rho_sand = ((2.6*0.7) + (1.0*0.3)).*rhow;
  Cp_rock = 0.2.*Cpw;
  Cp_sand = ((0.3*0.7) + (1.0*0.3)).*Cpw;

  rhob = (frac*rho_sand) + ((1-frac)*rho_rock);
  Cpb = (frac*Cp_sand) + ((1-frac)*Cp_rock);


  % Boundary layer depth - depth of substrate in thermal equilibrium with water
  hbl = get_opt(opts,'benthic_boundary_layer_depth',[]);
  if ( isempty(hbl) )
    hbl = 0.03; %[m]
  end;
  if ( doDebug )
    disp([mfilename ': Using benthic boundary layer depth ' num2str(hbl*100) 'cm']);
  end;


  % Conversion factor for benthic temperature change: [W/m2] -> [K/hr]
  Kperhour = 3600./(rhob*Cpb*hb);


  %% Parameters for benthic long-wave heat flux
  sigma = 5.67e-8;  % Stefan-Boltzman constant			%[W/m^2/K^4]

  % Seawater emissivity (to match STATION_BULK_ULR: see Anderson, 1952)
  % Xue et al (1998) say: "0.97 after Anderson (1952) and Reed (1976)"
  % Kraus and Businger suggest 0.98
  %epsw = 0.97;
  epsw = get_opt(opts,'b_epsw',0.97);

  % Benthos emissivity
  % % Taken from Evans et al. (1998), cites Oke (1978).
  % epsb = 0.98;
  % Sand~0.76, limestone~0.95 - saturated substrate may be closer to 0.89?
  % Also as above assume marine sediment=70% sandy clay+30% water: eps~0.82
  eps_sand = ((0.76*0.7) + (epsw*0.3));
  eps_other = 0.98;
  epsb = (frac*eps_sand) + ((1-frac)*eps_other);


  %% Parameters for benthos-water convective heat flux
  % Turbulent heat transfer coefficient: Cd^2 per Thibodeaux&Boyle (1987)
  % or variously Cd~2.5e-3 (Gross & Nowell 1983); 0.9-1.5e-3 (the range in
  % Reidenbach et al 2006); or presumed 0.017^2 (from Davis&Monismith 2011).
  Cbd = get_opt(opts,'b_convective_coefficient',[]);
  if ( isempty(Cbd) )
    %Cbd = 3.8e-4;
    %Cbd = (0.017*0.017);
    %Cbd = (0.017*1e-3);
    Cbd = (0.017*1e-2);
  end;
  if ( doDebug )
    disp([mfilename ': Using convective coefficient Cd^2~' num2str(Cbd)]);
  end;

  % NOTE: Tidal current may not dominate mean current: LONF1 TMD tide speed
  % has R^2~0.99, RMSE~0.06, slope 1 bias 0 vs. combined tide + FKEYS HYCOM.
  % However, at MLRF, R^2 drops to 0.27, RMSE~0.35, slope 0.37 bias 0.44.
  spd = uv_to_spd(u,v);	% Tidal/mean current above log layer	%[m/s]


  % Parameters for benthic heat conduction

  % Benthic thermal conductivity, mass avg.: marine sediment (1; Nobes et al
  % 1986) and porous water-sorbed carbonate (2.8; Thomas, Frost, Harvey 1973)
  Kb_sand = 1;
  Kb_rock = 2.8;
  Kb = (frac*Kb_sand) + ((1-frac)*Kb_rock);			%[W/m/K]


  % Temperature at base of thermally active sediment/rock layer: based on
  % annual mean sea temperatures at LONF1,MLRF1,SMKF1,FWYF1 for 1992-2010,
  % according to the method of Golosov and Kirillin (2010)

  % % General mean over all SEAKEYS sites
  % tbase = 26.5;							%[oC]

  tbase = get_opt(opts,'thermal_layer_base_temp',[]);
  if ( isempty(tbase) )
    tbase = nanmean(t);							%[oC]
    if ( doDebug )
      disp([mfilename ': Using thermal layer base temp T_base~Mean(T_s)']);
    end;
  end;


  % Benthos temperature
  tb = repmat(nan,size(t));
  %DEBUG:  disp(size(tb));

  % Net heat flux from benthos into the water column
  qbo = repmat(nan,size(t));


  % "Endo-Upwelling" aka heat flux through convective circulation in porous
  % substrate: According to estimates in Jean-Baptiste and Leclerc (2000),
  % this is unlikely to be a significant term in the reef heat budget.


  % We have no sea-bed temperature data, so we assume an initial temperature
  % equilibrium between the sea-floor and the water column

  % Assume only an initial equilibrium (free-running heat-exchange model)
  gapix = [1 ; length(dts)+1];

  % % Assume initial equilibrium at the start of any >1 hour gap in data
  % gapix = [1 ; (find(diff(dts) > (1.1/24))+1) ; length(dts)+1];

  % % Assume initial equilibrium at the start of any multi-day gap in data
  % gapix = [1 ; (find(diff(dts) > 1)+1) ; length(dts)+1];

  % % Assume on January 01 each year, sea and benthos at equilibrium
  % gapix = [1 ; find(get_yearday(dts)==0 & get_hour(dts)==0) ; length(dts)+1];

  % % Assume at local dawn each day [=.10 GMT], sea and benthos at equilibrium
  % gapix = [1 ; find(get_hour(dts)==10) ; length(dts)+1];

  % Pre-allocate arrays for speed
  qbswi = repmat(nan,size(tb));
  qblwi = repmat(nan,size(tb));
  qblwo = repmat(nan,size(tb));
  qbsh = repmat(nan,size(tb));
  qbcdu = repmat(nan,size(tb));
  qbcdd = repmat(nan,size(tb));

  %DEBUG:  disp(mfilename); disp(length(gapix)-1);
  for gapixix = 1:length(gapix)-1
    contigix = gapix(gapixix):(gapix(gapixix+1)-1);
    tb(contigix(1)) = t(contigix(1));

    % Benthic absorbed insolation - calculated as a by-product based on
    % sea-surface insolation, attenuation, and bottom reflectance in, e.g.,
    % STATION_ABSORBED_INSOLATION (v.)
    qbswi(contigix) = qb(contigix);

    qblwi(contigix) = epsw.*sigma.*( (t(contigix)+273.14) .^4 );

    %DEBUG:    disp(length(contigix));
    %DEBUG:    nby10 = floor(length(contigix)/10);
    for ix = 1:length(contigix)
      % Estimate outgoing long-wave flux for each hour
      qblwo(contigix(ix)) = epsb.*sigma.*( (tb(contigix(ix))+273.14) .^4 );
      % Estimate convective heat flux from benthos into water column
      %%%%???DEBUG Oops! Do we multiply first by rhow*Cpw, then by rhob*Cpb*hb below?
      qbsh(contigix(ix)) = rhow.*Cpw.*Cbd.*spd(contigix(ix)).*(t(contigix(ix)) - tb(contigix(ix)));

      % Estimate heat conduction across sea-bed/water interface upward into water
      qbcdu(contigix(ix)) = -Kb * ((t(contigix(ix)) - tb(contigix(ix)))/hbl);

      % Estimate heat conduction through sediment/coral rock layer downward
      %   Qcd = K(dT/dz), see definitions of KB,TBASE above for assumptions
      qbcdd(contigix(ix)) = -Kb * ((tb(contigix(ix)) - tbase)/hb);

      % Estimate next hour's benthic temperature based on total fluxes
      tb(contigix(ix)+1) = tb(contigix(ix)) + ...
          (Kperhour.*( qbswi(ix) - qblwo(ix) + qblwi(ix) + qbsh(ix) - qbcdu(ix) + qbcdd(ix) ));
      %DEBUG:      if (mod(ix,nby10)==0); disp(ix); end;
    end;

    % Estimate heat flux INTO water column
    qbo(contigix) = - qblwo(contigix) + qblwi(contigix) + qbsh(contigix) - qbcdu(contigix);

    % ??? DEBUG:
    % At each transition (gap, dawn, year-end), rebalance energy budget with
    % an UNPHYSICAL radiative "pulse" from the benthos into the water column!
    if ( contigix(end) < length(dts) )
      qbo(contigix(end)) = (tb(contigix(end))-t(contigix(end)+1))./Kperhour;
      if ( doDebug )
        disp([mfilename ': RESET T_B-T_s=0 at ' datestr(dts(contigix(end)))]);
      end;
    end;
  end;

  tb(length(dts)) = [];

  stn.(tbfld).date = dts;
  stn.(tbfld).data = tb;
  stn.(qbofld).date = dts;
  stn.(qbofld).data = qbo;

  % DEBUG: If caller asks for diagnostic time series, return them
  if ( exist('diagfld','var') && ~isempty(diagfld) )
    stn.(diagfld).qbswi.date = dts(:);
    stn.(diagfld).qbswi.data = qbswi(:);
    stn.(diagfld).qblwi.date = dts(:);
    stn.(diagfld).qblwi.data = qblwi(:);
    stn.(diagfld).qblwo.date = dts(:);
    stn.(diagfld).qblwo.data = qblwo(:);
    stn.(diagfld).qblw.date = dts(:);
    stn.(diagfld).qblw.data = qblwi(:) - qblwo(:);
    stn.(diagfld).qbsh.date = dts(:);
    stn.(diagfld).qbsh.data = qbsh(:);
    stn.(diagfld).qbcdd.date = dts(:);
    stn.(diagfld).qbcdd.data = qbcdd(:);
    stn.(diagfld).qbcdu.date = dts(:);
    stn.(diagfld).qbcdu.data = qbcdu(:);
    stn.(diagfld).qbcd.date = dts(:);
    stn.(diagfld).qbcd.data = qbcdd(:) + qbcdu(:);
  end;

  %DEBUG:  fmg; plot(stn.(diagfld).qbsh.date,-cumsum(stn.(diagfld).qbsh.data-stn.(diagfld).qbcdu.data+stn.(diagfld).qblwi.data-stn.(diagfld).qblwo.data),'r'); plot(stn.(diagfld).qbswi.date,cumsum(stn.(diagfld).qbswi.data),'b'); legend('Q_b_S_H+Q_b_C_D^U+Q_b_L_W','Q_b_S_W^I'); datetick3;

  %DEBUG:  toc,

return;
