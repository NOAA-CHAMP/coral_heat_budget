function res = horizontal_convection(t,s,h,q0,bet,opts,dts,q0lp,wndlp)
%function res = horizontal_convection(t,s,h,q0,bet,opts,dts,q0lp,wndlp)
%
% Calculate "horizontal convective" (Monismith et al 2006) velocity and heat
% advection, based on vectors of sea temperature (T), salinity (S), site
% depth [m] (H), surface air-sea flux [W/m^2] (Q0), and site bathymetryic
% slope or 'beta' (BET). Also accepts optional struct OPTS (see GET_OPT).
% Returns a struct of results, RES, with fields for: buoyancy flux .B0;
% convective velocity scale .uf; volumetric flow .Qv_* and velocity .u_*
% under four Monismith et al scalings (2006); estimated true convective
% velocity .u; differential heating terms .dTd*; convective exchange rate
% .dTdthc; and total heat storage change per hour .dTdt. 
%
% Last Saved Time-stamp: <Sun 2014-09-28 18:30:18 Eastern Daylight Time gramer>

  doDebug = get_opt(opts,'hc_debug',true);

  rho = sw_dens(s,t,h);						%[kg/m^3]
  Cp = sw_cp(s,t,h);						%[J/kg*K]
  rhoCp = rho.*Cp;						%[J/K*m^3]
  g = 9.79;							%[m/s^2]
  alph = sw_alpha(s,t,h);					%[1/K]

  dt = 3600;							%[s/hr]

  fac = dt./(rhoCp.*h);						%[K*m^2/W*hr]
  res.termFactor = fac;

  % Buoyancy flux
  res.q0 = q0;
  res.B0=(g.*alph.*abs(q0))./(rhoCp);				%[m^2/s3]

  % Characteristic convective velocify u_f:
  %  Fischer et al (1979), Eq. 6.40 - per Monismith et al (2006)
  res.uf=(res.B0.*h).^(1/3);					%[m/s]


  % Periodicity of thermal forcing (DEFAULT: diurnal warming, 24 h)
  Tf = get_opt(opts,'hc_Tf',24);

  % Volumetric flow: Panel letter refers to Fig. 10, Monismith et al (2006)

  % Advective inertial momentum balance, steady thermal balance: Panel (c)
  res.Qv_SS=res.uf.*h./(bet.^(1/3));				%[m^2/s]

  % Viscous/unsteady inertia, balanced thermal forcing: Panel (a)
  %res.Qv_US=sqrt((res.uf.^3) .* (Tf.*dt) ./ h);		%[m/s???]
  res.Qv_US=sqrt((res.uf.^3) .* (Tf.*dt) .* h);			%[m^2/s]

  % Advective inertial balance, unbalanced thermal forcing: Panel (f)
  %res.Qv_SU=(res.uf.*(Tf.*dt)./h).^(3/2);			%[1???]
  % VNU -> VNS: (uf^3*T/D)^1/2 * bet^-1/3 * (D/T*uf)^1/2 = uf*bet^-1/3
  % VUU -> VUS: bet*uf^3(T^2/D^2) * bet^-1/3 * (D/T*uf)^1/2 = bet^2/3*uf*(uf*T/D)^3/2
  res.Qv_SU=(bet.^(2/3)).*res.uf.*((res.uf.*(Tf.*dt)./h).^(3/2));		%[m^2/s]

  % Viscous/unsteady inertia, unbalanced thermal forcing: Panel (d)
  %res.Qv_UU=bet.*(res.uf.^3).*((Tf.*dt)^2)./(h.^2);		%[m/s???]
  res.Qv_UU=bet.*(res.uf.^3).*((Tf.*dt)^2)./h;			%[m^2/s]


  % Assuming fit in Fig. 10, panel (c) of Monismith et al (2006)
  res.u_SS=(5.0*res.Qv_SS./h) - 0.05;   res.u_SS(res.u_SS<0) = 0;
  % Assuming fit in Fig. 10, panel (a) of Monismith et al (2006)
  res.u_US=(3.0*res.Qv_US./h) - 0.026;  res.u_US(res.u_US<0) = 0;
  % Assuming fit in Fig. 10, panel (f) of Monismith et al (2006)
  res.u_SU=(2.7*res.Qv_SU./h) - 0.0224; res.u_SU(res.u_SU<0) = 0;
  % Assuming fit in Fig. 10, panel (d) of Monismith et al (2006)
  res.u_UU=(0.1*res.Qv_UU./h);	        res.u_UU(res.u_UU<0) = 0;

  % Introduce convective response lags to (unsteady) periodic forcing
  lagoff = get_opt(opts,'hc_convective_lag',3);
  res.u_SU(lagoff+1:end) = res.u_SU(1:end-lagoff);  res.u_SU(1:lagoff) = 0;
  res.u_UU(lagoff+1:end) = res.u_UU(1:end-lagoff);  res.u_UU(1:lagoff) = 0;

  % Net transport ex mixing efficiency (DEFAULT: 8% mixing efficiency)
  R = get_opt(opts,'hc_R',(1.00-0.08));

  sclg = get_opt(opts,'hc_scaling','SS');

  descStr = [upper(mfilename),' ',sclg];

  res.scaling = sclg;
  switch ( upper(sclg) ),
   case 'SS',  res.u = res.u_SS;
   case 'US',  res.u = res.u_US;
   case 'SU',  res.u = res.u_SU;
   case 'UU',  res.u = res.u_UU;

   case 'SEASONAL_ORIGINAL_RECIPE',
    % Per-month scaling:
    %    1-3: US, R=1.00
    %    4-6: UU, R=0.92
    %    7-9: SU, R=1.00
    %  10-12: SS, R=0.92
    res.u = res.u_SS;
    R = repmat(R,size(res.u));
    seasix = find(get_season(dts)==1);
    res.u(seasix) = res.u_US(seasix);	R(seasix) = 1.00;
    seasix = find(get_season(dts)==2);
    res.u(seasix) = res.u_UU(seasix);	R(seasix) = 0.92;
    seasix = find(get_season(dts)==3);
    res.u(seasix) = res.u_SU(seasix);	R(seasix) = 1.00;
    seasix = find(get_season(dts)==4);
    res.u(seasix) = res.u_SS(seasix);	R(seasix) = 0.92;

   case 'STEADY_SEASONAL',
    res.u = res.u_SS;
    R = repmat(R,size(res.u));
    seasix = find(ismember(get_month(dts),[1:4,11:12]));
    res.u(seasix) = res.u_SU(seasix);

   case 'UNSTEADY_SEASONAL',
    res.u = res.u_US;
    R = repmat(R,size(res.u));
    seasix = find(ismember(get_month(dts),[1:4,11:12]));
    res.u(seasix) = res.u_UU(seasix);


   case 'SEASONAL_STEADY',
    res.u = res.u_SS;
    R = repmat(R,size(res.u));
    seasix = find(ismember(get_month(dts),[1:4,11:12]));
    res.u(seasix) = res.u_US(seasix);

   case 'SEASONAL_UNSTEADY',
    res.u = res.u_SU;
    R = repmat(R,size(res.u));
    seasix = find(ismember(get_month(dts),[1:4,11:12]));
    res.u(seasix) = res.u_UU(seasix);


   case 'SEASONAL_SEASONAL',
    res.u = res.u_US;
    R = repmat(R,size(res.u));
    seasix = find(ismember(get_month(dts),[1:4,11:12]));
    res.u(seasix) = res.u_SU(seasix);


   case 'ADAPTIVE_ORIGINAL_RECIPE',
    res.u = res.u_US;
    R = repmat(R,size(res.u));

    % typix = find(wndlp  < 5 & abs(q0lp)  < 300);
    % res.u(typix) = res.u_SS(typix);	R(typix) = 1.00;

    typix = find(wndlp >= 5 & abs(q0lp)  < 300);
    res.u(typix) = res.u_US(typix);	R(typix) = 1.00;

    % typix = find(wndlp  < 5 & abs(q0lp) >= 300);
    % res.u(typix) = res.u_SU(typix);	R(typix) = 0.92;

    typix = find(wndlp >= 5 & abs(q0lp) >= 300);
    res.u(typix) = res.u_UU(typix);	R(typix) = 0.92;

   case 'ADAPTIVE',
    res.u = res.u_SU;

    typix = find(wndlp >= 10);
    res.u(typix) = res.u_UU(typix);
  end;

  % Reality check
  res.u(res.u<0) = 0;

  % Adjust algorithm based on cooling/heating and other factors
  coolt = get_opt(opts,'hc_cooling_cutoff',0);
  warmt = get_opt(opts,'hc_warming_cutoff',0);
  coolix = find(q0<coolt);
  zeroix = find(coolt<=q0 & q0<=warmt);
  warmix = find(q0>warmt);

  res.hc_cooling_factor = get_opt(opts,'hc_cooling_factor',1.0);
  res.hc_warming_factor = get_opt(opts,'hc_warming_factor',1.0);
  res.u(coolix) = res.hc_cooling_factor*res.u(coolix);
  res.u(zeroix) = 0;
  res.u(warmix) = res.hc_warming_factor*res.u(warmix);


  % Onset delay for convective circulation
  % Chubarenko 2010
  l = dt .* res.uf;
  D = h + (bet.*l);
  dT = dt .* ( (q0./(rhoCp.*(h+D))) - (q0./(rhoCp.*h)) );
  % dRho = dT .* alpha,
  rhox = sw_dens(s,(t+dT),h);
  dRho = (rho - rhox) ./ rhox;
  res.onT = sqrt(l./(g.*abs(dRho).*abs(bet)));


  % Remove any unphysical exchange velocity estimates
  maxu = get_opt(opts,'hc_maximum_velocity',2.00);
  badix = find(res.u>maxu);
  if (doDebug)
    disp([descStr ': Zeroing ' num2str(length(badix)) ' unphysical speeds >' num2str(maxu) 'm/s']);
  end;
  res.u(badix) = 0;
  if (doDebug)
    %disp([descStr ': Truncating ' num2str(length(badix)) ' unphysical speeds >' num2str(maxu) 'm/s']);
  end;
  % res.u(badix) = maxu;

  % If onset time is too long, assume thermal siphon simply does not operate
  maxt = get_opt(opts,'hc_max_onset_secs',dt);
  badix = find(res.onT>=maxt);
  if (doDebug)
    disp([descStr ': Zeroing ' num2str(length(badix)) ' convective onsets >' num2str(maxt) 's']);
  end;
  res.u(badix) = 0;

  % Calculate resultant heat advection
  res.dx=res.u.*dt;
  res.dTdtq0=dt.*q0./(rhoCp.*h);
  res.dTdtx=dt.*q0./(rhoCp.*(h+(bet.*res.dx)));
  res.dTdx=(res.dTdtq0-res.dTdtx)./res.dx;
  res.dTdx(~isfinite(res.dTdx)) = 0;
  res.dTdthc=-R.*dt.*res.u.*res.dTdx;

  res.dTdt=res.dTdtq0+res.dTdthc;


return;
