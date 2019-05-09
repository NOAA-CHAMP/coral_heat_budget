function stn = station_heat_budget_errors(stn,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,substitute_field_names)
%function stn = station_heat_budget_errors(stn,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,substitute_field_names)
%
% Calculate errors for each term in reef ocean heat budget. Uses TOGA/COARE
% error algorithm for turbulent fluxes, and published relative and absolute
% errors, confirmed by local in situ regressions where ever available, to
% calculate error budget for all other terms. See STATION_HEAT_BUDGET or
% OPTIM_Q0 for args. Adds new fields to STN with '_err' appended to name.
%
% Last Saved Time-stamp: <Fri 2012-08-03 16:59:40  lew.gramer>


  set_more off
  %tic,

  if ( ~isfield(stn,'station_name') )
    error('Station struct STN had no .station_name field!');
  end;

  sal = get_opt(stn.opts,'default_salinity',35.5);


  %%%
  %% Call SCRIPT to set:
  %% Set variable-name prefixes ("PFX") for various input/output datasets;
  %% AND, set all station struct fieldnames used to produce heat budget 
  %% NOTE WELL: Calls FIX_VARNAMELENGTHS to fix variable name string lengths
  %%             to meet MATLAB limitations. MAJOR SIDE EFFECTS...
  station_heat_budget_field_names;


  % Needed for TOGA/COARE error algorithms
  if ( isempty(strfind(path,'fairall')) )
    FAIRALLHOME = get_ecoforecasts_path('../fairall');
    addpath(FAIRALLHOME);
    clear FAIRALLHOME
    rehash
  end;


  [zu, zt] = station_instrument_heights(stn.station_name);
  zq = zt;


  [dix,six,aix,qaix,qsix,Wix,hix,lix,Hsix,qlhix,qshix,qswiix,qswix,qlwiix,...
   sq0ix,sqtix,qbix,qbtix,bq0ix,bq0tix,fqudTix,kd2Tix,hchcix,hcix] = ...
      intersect_all_dates([],...
                          stn.(diagfld).date,...
                          stn.(sfld).date,...
                          stn.(afld).date,...
                          stn.(qafld).date,...
                          stn.(qsfld).date,...
                          stn.(Wfld).date,...
                          stn.(mhfld).date,...
                          stn.(pblzfld).date,...
                          stn.(whfld).date,...
                          stn.(qlhfld).date,...
                          stn.(qshfld).date,...
                          stn.(dsrfld).date,...
                          stn.(srfld).date,...
                          stn.(dlrfld).date,...
                          stn.(sq0fld).date,...
                          stn.(sqtfld).date,...
                          stn.(qbofld).date,...
                          stn.(qbotfld).date,...
                          stn.(bq0fld).date,...
                          stn.(bq0tfld).date,...
                          stn.(fqudTfld).date,...
                          stn.(kd2Tfld).date,...
                          stn.(hcdTdthc).date,...
                          stn.(hcdTdt).date);

  dts = stn.(sfld).date(six);

  s = stn.(sfld).data(six);
  a = stn.(afld).data(aix);
  qss = stn.(qsfld).data(qsix);
  qsa = stn.(qafld).data(qaix);
  W = kts2mps(stn.(Wfld).data(Wix));
  h = stn.(mhfld).data(hix);
  pblz = stn.(pblzfld).data(lix);
  Hs = stn.(whfld).data(Hsix);
  qlh = real(stn.(qlhfld).data(qlhix));	qlh(~isfinite(qlh)) = 0;	qlh2 = signedSQR(qlh);
  qsh = real(stn.(qshfld).data(qshix));	qsh(~isfinite(qsh)) = 0;	qsh2 = signedSQR(qsh);
  qswi = stn.(dsrfld).data(qswiix);	qswi(~isfinite(qswi)) = 0;	qswi2 = signedSQR(qswi);
  qsw = stn.(srfld).data(qswix);	qsw(~isfinite(qsw)) = 0;	qsw2 = signedSQR(qsw);
  qlwi = stn.(dlrfld).data(qlwiix);	qlwi(~isfinite(qlwi)) = 0;	qlwi2 = signedSQR(qlwi);
  sq0 = stn.(sq0fld).data(sq0ix);  	sq0(~isfinite(sq0)) = 0;  	sq02 = signedSQR(sq0);
  sqt = stn.(sqtfld).data(sqtix);  	sqt(~isfinite(sqt)) = 0;  	sqt2 = signedSQR(sqt);
  qbo = stn.(qbofld).data(qbix);  	qbo(~isfinite(qbo)) = 0;  	qbo2 = signedSQR(qbo);
  qbot = stn.(qbotfld).data(qbtix);  	qbot(~isfinite(qbot)) = 0;  	qbot2 = signedSQR(qbot);
  bq0 = stn.(bq0fld).data(bq0ix);	bq0(~isfinite(bq0)) = 0;	bq02 = signedSQR(bq0);
  bq0t = stn.(bq0tfld).data(bq0tix);	bq0t(~isfinite(bq0t)) = 0;	bq0t2 = signedSQR(bq0t);
  fqudT = stn.(fqudTfld).data(fqudTix);	fqudT(~isfinite(fqudT)) = 0;	fqudT2 = signedSQR(fqudT);
  kd2T = stn.(kd2Tfld).data(kd2Tix);	kd2T(~isfinite(kd2T)) = 0;	kd2T2 = signedSQR(kd2T);
  hchc = stn.(hcdTdthc).data(hchcix);	hchc(~isfinite(hchc)) = 0;	hchc2 = signedSQR(hchc);
  hc = stn.(hcdTdt).data(hcix); 	hc(~isfinite(hc)) = 0;  	hc2 = signedSQR(hc);

  switch ( TIDEPFX ),
   case 'tmd_tide',
    % TMD Tide error based on SMKF1 tide sensor: regression slope 0.7, RMSE 0.1[m]
    %  SMKF1 bias was 0.75[m]; also TRNF1 slope 1.01, RMSE 0.1[m], bias -1.1[m]
    % This is needed both for Absorbed Insolation (gamma*Qsw) and Stored Heat errors
    s2h = 0.1 + (0.30.*h);					sh = signedsqrt(s2h);
   case 'tpxo_tide',
    % TPXO Tide error based on SMKF1 tide sensor:
    s2h = 0.1 + (0.30.*h);					sh = signedsqrt(s2h);
   otherwise,
    error('Unknown tide source "%s"',TIDEPFX);
  end;  


  switch ( WAVEPFX ),
   % NOTE: Site 4 in 8.8m on crest - other two are in deeper water offshore
   case 'erai',
    % ERAI significant wave height error AFTER empirical corrections using
    % Haus study sites: Site 3 slope 1.04, bias 0.01, r^2=0.46, RMSE~0.21;
    % Site 4: 0.87, -0.08, 0.35, 0.23; Site 7: 0.97, -0.00, 0.46, 0.21.
    s2Hs = 0.21 + (0.10.*Hs);					sHs = signedsqrt(s2Hs);
   case 'ww3',
    % WW3 significant wave height error AFTER empirical corrections using
    % Haus study sites: Site 3 slope 1.01, bias 0.01, r^2=0.72, RMSE~0.15;
    % Site 4: 0.84, -0.07, 0.44, 0.22; Site 7: 1.03, -0.05, 0.74, 0.15.
    s2Hs = 0.20 + (0.15.*Hs);					sHs = signedsqrt(s2Hs);
   case 'ndbc',
    % Error for bulk significant wave height from NDBC wind data using
    % Haus study sites: Site 3 slope 1.00, bias 0.00, r^2=0.50, RMSE~0.20;
    % Site 4: 0.74, -0.07, 0.10, 0.24; Site 7: 1.03, -0.07, 0.41, 0.22.
    % error('NEED TO VERIFY ERROR SLOPE FOR WAVEPFX="ndbc"!');
    % s2Hs = 0.23 + (0.22.*Hs);					sHs = signedsqrt(s2Hs);
    s2Hs = 0.25 + (0.35.*Hs);					sHs = signedsqrt(s2Hs);
   otherwise,
    error('Unknown wave source "%s"',WAVEPFX);
  end;  


  % Adjust heat storage errors for uncertainty of T (~0.2K), S (~1PSU), and h (~3m) in K/hr calculation
  % h=10; t=25; s=35.5;  sh=3; st=0.2; ss=1.0;  hp=h+sh; tp=t+st; sp=s-ss;
  % sKperHr=3600.*( (1/(h.*sw_cp(s,t,h).*sw_dens(t,s,h)))-(1/(hp.*sw_cp(sp,tp,hp).*sw_dens(tp,sp,hp))) );
  sKperHr = 1e-5;
  sWm2 = (1/sKperHr);


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% RADIATIVE FLUXES

  % Estimate shortwave radiative flux errors

% ??? VERIFY against Lewis, Wei, van Dommelen, Voss (2011)
% ??? (see same comment in STATION_ABSORBED_INSOLATION.m)

  % % Consistent with Markovic et al (2009), 5% relative error, 7 W/m^2 bias in ERA40 Qswi
  % sQswi = 7 + (0.05.*qswi);					s2Qswi = signedSQR(sQswi);

  % % Based on regressing FWYF1 ERAI vs. Lew QC'd RSMAS in situ DAILY AVERAGES: bias -3, RMSE 40
  % sQswi = 43 + (0.01.*qswi);					s2Qswi = signedSQR(sQswi);

  % % Based on regressing ADJUSTED VKAF1 ERAI vs. Lew QC'd RSMAS in situ DAILY AVERAGES: bias +1, RMSE 35
  % sQswi = 35 + (0.004.*qswi);					s2Qswi = signedSQR(sQswi);

  % Mean of Markovic et al (2009) and regression: 2% relative error, RMSE 17, 3 W/m^2 bias
  sQswi = 20 + (0.02.*qswi);					s2Qswi = signedSQR(sQswi);

  % Remove "error" from night-time values
  sQswi(qswi<1) = 0;						s2Qswi = signedSQR(sQswi);


  % Assuming Albedo has constant relative uncertainty of 4% *of Qswi*
  % sQsw = 43 + (0.05.*qswi);					s2Qsw = signedSQR(sQsw);
  % sQsw = 36 + (0.044.*qswi);					s2Qsw = signedSQR(sQsw);
  % sQsw = 7 + (0.05.*qswi);					s2Qsw = signedSQR(sQsw);
  sQsw = 20 + (0.064.*qswi);					s2Qsw = signedSQR(sQsw);

  % Remove "error" from night-time values
  sQsw(qswi<1) = 0;						s2Qsw = signedSQR(sQsw);


  % Assuming gamma is time-varying and with constant 7% relative uncertainty
  % See ESTIM_GAMMA_ERROR.m for the basis of this relative error estimate
  [gamix,srfix] = intersect_dates(stn.(gamfld).date,stn.(srfld).date(qswix));
  gamma = stn.(gamfld).data(gamix);
  gamma(~isfinite(gamma)) = mean(gamma(isfinite(gamma)));
  gamma2 = signedSQR(gamma);

  sgamma = 0.07 .* gamma;					s2gamma = signedSQR(sgamma);

  [ig,r_sgamma_sQsw] = cov_ts(sgamma,sQsw,@isfinite,@isfinite);

  s2gammaQsw = ( (gamma2.*s2Qsw) + (qsw2.*s2gamma) + (2.*gamma.*qswi.*sgamma.*sQsw.*r_sgamma_sQsw)  );
  sgammaQsw = signedsqrt(s2gammaQsw);

  %DEBUG:
  s2gammaQsw1 = ( (gamma2.*s2Qsw)                                                                   );
  %DEBUG:
  s2gammaQsw2 = (                   (qsw2.*s2gamma)                                                 );
  %DEBUG:
  s2gammaQsw3 = (                                     (2.*gamma.*qswi.*sgamma.*sQsw.*r_sgamma_sQsw) );


  % Remove "error" from night-time values
  s2gammaQsw(qswi<1) = 0;					sgammaQsw = signedsqrt(s2gammaQsw);


  %%%% HACK??? (For Albedo and Gamma calculated from in situ data)
  %%%% DEBUG:  s2Qsw = signedSQR(1 - prctile(Albedo,99)) .* s2Qswi;
  %%%% DEBUG:  s2Qsw(qswi<0.1) = 0;
  %%%% DEBUG:  sQsw = signedsqrt(s2Qsw);
  %%%% DEBUG:  s2gammaQsw = signedSQR(nanmedian(gamma(0<gamma&gamma<1))) .* s2Qsw;	sgammaQsw = signedsqrt(s2gammaQsw);


  % Estimate longwave radiative flux errors

  % % Consistent with Markovic et al (2009), 1% relative error, 3 W/m^2 bias in ERA40 QlwI
  % sQlwi = 3 + (0.01.*qlwi);					s2Qlwi = signedSQR(sQlwi);
  % % Based on regressing FWYF1 ERAI vs. Lew QC'd RSMAS in situ: bias -10, RMSE 14
  % sQlwi = 24 + (0.021.*qlwi);					s2Qlwi = signedSQR(sQlwi);
  % % Based on regressing FWYF1 ERAI (after empirical adjustment) vs. Lew QC'd RSMAS in situ, R^2~0.84
  % sQlwi = 14 + (0.005.*qlwi);					s2Qlwi = signedSQR(sQlwi);

  % % Based on regressing ADJUSTED VKAF1 ERAI (after empirical adjustment) vs. Lew QC'd RSMAS in situ, R^2~0.93
  % sQlwi = 10 + (0.02.*qlwi);					s2Qlwi = signedSQR(sQlwi);

  % Mean of Markovic et al (2009) and regression: 1.5% relative error, RMSE 5, 1.5 W/m^2 bias
  sQlwi = 6.5 + (0.015.*qlwi);					s2Qlwi = signedSQR(sQlwi);

  % ULRF error est. based on graybody calculation with nominal 0.11K Ts error
  sQlwo = 1;							s2Qlwo = signedSQR(sQlwo);

  % Sum of variances (covariance between sQlwi and sQlwo ignored)
  s2Qlw = (s2Qlwi + s2Qlwo);					sQlw = signedsqrt(s2Qlw);


  % Save for future reference
  stn.([srfld '_err']).date = dts(:);
  stn.([srfld '_err']).data = sQsw(:);
  stn.([lrfld '_err']).date = dts(:);
  stn.([lrfld '_err']).data = sQlw(:);

  stn.([asrfld '_err']).date = dts(:);
  stn.([asrfld '_err']).data = sgammaQsw(:);


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% TURBULENT FLUXES

  %error specs  [intercept error   slope error]
  %Xtrue=+-interror+(1+-sloperror)*Xmeas
  ue=[.55 .03];%'true' wind speed relative to sea surface
  % .03 slope error verified by duplicate wind data at MLRF1,FWYF1,SMKF1
  te=[.09 0];%sea-air temperature difference
  qe=[0.18 0];%sea-air spechumid difference


  bet = 1.25;
  k = 0.4;
  al = 9.8 ./ (a + 273.15);

  qa = qsa .* 1e3;  % Algorithm expects [g/kg]
  qs = qss .* 1e3;  % Algorithm expects [g/kg]

  zi = pblz;
  dq = qs - qa;
  dt = s - a;
  dth = dt + ( 6.1e-4 .* (a + 273.15) .* dq );

  eu = sqrt( (ue(1).^2) + ((ue(2).*W).^2) );
  et = sqrt( (te(1).^2) + ((te(2).*dt).^2) );
  eq = sqrt( (qe(1).^2) + ((qe(2).*dq).^2) );

  ecu = 0.001 ./ 10;
  ect = 0.001 ./ 10;
  ecq = 0.001 ./ 10;

  hs = stn.(qshfld).data(qshix);	hs(~isfinite(hs)) = 0;
  hl = stn.(qlhfld).data(qlhix);	hl(~isfinite(hl)) = 0;

  % COARE algorithms a bit sloppy - sometime produce complex numbers
  hs = real(hs);
  hl = real(hl);

  hb = hs + (6.1e-4 .* 300 .* hl ./ 2.5);
  L = stn.(diagfld).monin(dix);
  zetu = zu ./ L;
  zett = zt ./ L;
  zetq = zq ./ L;

  ztn = 10 ./ L;
  if ( isfield(stn.(diagfld),'z0u') )
    cunh = k ./ log(zu ./ stn.(diagfld).z0u(dix));
    ctnh = k ./ log(zt ./ stn.(diagfld).z0t(dix));
    cqnh = k ./ log(zq ./ stn.(diagfld).z0q(dix));
  else
    z0u = zu;  %%%%???DEBUG: Probably not valid!
    cunh = k ./ log(zu ./ z0u);
    ctnh = k ./ log(zt ./ z0u);
    cqnh = k ./ log(zq ./ z0u);
  end;

  dcu = (psiu_30(zetu) - psiu_30(0.99.*zetu)) ./ 0.01;
  dct = (psit_30(zett) - psit_30(0.99.*zett)) ./ 0.01;
  dcq = (psit_30(zetq) - psit_30(0.99.*zetq)) ./ 0.01;

  usr = stn.(diagfld).ustar(dix);
  tsr = stn.(diagfld).tstar(dix);
  qsr = stn.(diagfld).qstar(dix);

  ws = ( al .* usr .* abs(tsr + (6.1e-4 .* 300 .* qsr)) .* zi ) .^ 0.333;

  su = (1 - (cunh ./ k .* psiu_30(zetu)));
  st = (1 - (ctnh ./ k .* psit_30(zett)));
  sq = (1 - (cqnh ./ k .* psit_30(zetq)));
  uw = sqrt( (W.*W) + ((bet.*ws).^2) );
  aw = ( (bet.*ws) ./ uw ) .^2;

  fs = cunh ./ k .* dcu ./ su + aw ./ 3;
  fu = cunh ./ k .* dcu ./ su;
  ft = ctnh ./ k .* dct ./ st;
  fq = cqnh ./ k .* dcq ./ sq;

  %***** Transfer coefficient error contribution
  xs1 = ecu ./ su ./ ctnh;
  xt1 = ect ./ st ./ ctnh;
  xq1 = ecq ./ sq ./ cqnh;
  xu1 = ecu ./ su ./ cunh;

  %***** Sensor error contribution
  xs2 = eu ./ W .* (W ./ uw) .^ 2;
  xt2 = et ./ dt;
  xq2 = eq ./ dq;
  xu2 = eu ./ W;

  %************** Analytical error functions
  d = (1 - ft) .* (1 - aw) + 2 .* fs; % Universal denominator
  dus = sqrt( (1 - ft) .^ 2 .* (xs1 .^ 2 + xs2 .^ 2) + fs .^ 2 .* (xt1 .^ 2 + xt2 .^ 2) ) ./ d;
  dzet = sqrt( (1 - aw) .^ 2 .* (xt1 .^ 2 + xt2 .^ 2) + 4 .* (xs1 .^ 2 + xs2 .^ 2) ) ./ d;


  sQsh = sqrt( ((1 - aw) + 3 .* fs) .^ 2 .* (xt1 .^ 2 + xt2 .^ 2) + (xs1 .^ 2 + xs2 .^ 2) .* (1 - 3 .* ft) .^ 2 ) ./ d;
  sQlh = sqrt( xq1 .^ 2 + xq2 .^ 2 + ((fq .* (1 - aw) + fs) .^ 2 .* (xt1 .^ 2 + xt2 .^ 2) + ((1 - ft) - 2 .* fq) .^ 2 .* (xs1 .^ 2 + xs2 .^ 2)) ./ d .^ 2 );
  sTauMaybe = sqrt( (xu1 .* ((1 - ft) .* (2 - aw) + .67 .* aw)) .^ 2 + (xu2 .* ((1 - ft) .* (1 + (W ./ uw) .^ 2 - aw) + 2 .* fu .* (1 - (W ./ uw) .^ 2) + .67 .* aw)) .^ 2 + (xt1 .^ 2 + xt2 .^ 2) .* (fu .* (2 - aw) + .33 .* aw) .^ 2 ) ./ d;

  % COARE algorithms a bit sloppy - sometime produce complex numbers
  sQsh = real(sQsh);
  sQlh = real(sQlh);

  % Gross quality control on the relative errors
  sQsh(~isfinite(sQsh)|0>sQsh|sQsh>1) = 1;
  sQlh(~isfinite(sQlh)|0>sQlh|sQlh>1) = 1;

  stn.([qshfld '_relerr']).date = dts(:);
  stn.([qshfld '_relerr']).data = sQsh(:);
  stn.([qlhfld '_relerr']).date = dts(:);
  stn.([qlhfld '_relerr']).data = sQlh(:);

  % Above calculates relative (normalized) errors - convert to absolute errors!
  sQsh = sQsh(:) .* abs(hs(:));					s2Qsh = signedSQR(sQsh);
  sQlh = sQlh(:) .* abs(hl(:));					s2Qlh = signedSQR(sQlh);


  stn.([qshfld '_err']).date = dts(:);
  stn.([qshfld '_err']).data = sQsh(:);
  stn.([qlhfld '_err']).date = dts(:);
  stn.([qlhfld '_err']).data = sQlh(:);


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% NET AIR-SEA FLUXES

  % Estimate error covariances for total Q0 error
  COVgammaQswQlw = cov_ts(sgammaQsw,sQlw,@isfinite,@isfinite);
  COVgammaQswQsh = cov_ts(sgammaQsw,sQsh,@isfinite,@isfinite);
  COVgammaQswQlh = cov_ts(sgammaQsw,sQlh,@isfinite,@isfinite);
  COVQswQlw = cov_ts(sQsw,sQlw,@isfinite,@isfinite);
  COVQswQsh = cov_ts(sQsw,sQsh,@isfinite,@isfinite);
  COVQswQlh = cov_ts(sQsw,sQlh,@isfinite,@isfinite);
  COVQlwQsh = cov_ts(sQlw,sQsh,@isfinite,@isfinite);
  COVQlwQlh = cov_ts(sQlw,sQlh,@isfinite,@isfinite);
  COVQlhQsh = cov_ts(sQlh,sQsh,@isfinite,@isfinite);

  s2Q0component = s2Qsw + s2Qlw + s2Qlh + s2Qsh;
  s2Q0component(~isfinite(s2Q0component)) = 0;			sQ0component = signedsqrt(s2Q0component);
  COVQlhQ0component = cov_ts(sQlh,sQ0component,@isfinite,@isfinite);

  s2Q0 = ( s2Q0component + ...
           (2.*(COVQswQlw + COVQswQlh + COVQswQsh + COVQlwQlh + COVQlwQsh + COVQlhQsh)) );
  s2Q0(~isfinite(s2Q0)) = 0;					sQ0 = signedsqrt(s2Q0);

  s2gQ0 = ( s2gammaQsw + s2Qlw + s2Qlh + s2Qsh + ...
            (2.*(COVgammaQswQlw + COVgammaQswQlh + COVgammaQswQsh + COVQlwQlh + COVQlwQsh + COVQlhQsh)) );
  s2gQ0(~isfinite(s2gQ0)) = 0;					sgQ0 = signedsqrt(s2gQ0);

  stn.([sq0fld '_err']).date = dts(:);
  stn.([sq0fld '_err']).data = sQ0(:);

  stn = station_heat_flux_term(stn,[sq0fld,'_err'],[sqtfld,'_err'],sfld,sal,mhfld);
  % Adjust for uncertainty of T (0.5K), S (1.0PSU), and h (18m) in K/hr calculation
  sQt = stn.([sqtfld,'_err']).data + (sKperHr.*abs(sq0));	s2Qt = signedSQR(sQt);
  stn.([sqtfld,'_err']).data(:) = sQt(:);

  stn.([sq0fld,'_qlh_coverr']).date = dts(:);
  stn.([sq0fld,'_qlh_coverr']).data = COVQlhQ0component(:);

  stn.([q0fld '_err']).date = dts(:);
  stn.([q0fld '_err']).data = sgQ0(:);


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% BENTHIC FLUX

  sQb = 0.03.*abs(qbo);						s2Qb = signedSQR(sQb);

  stn.([qbofld '_err']).date = dts(:);
  stn.([qbofld '_err']).data = sQb(:);

  stn = station_heat_flux_term(stn,[qbofld,'_err'],[qbotfld,'_err'],sfld,sal,mhfld);
  % Adjust for uncertainty of T (0.5K), S (1.0PSU), and h (18m) in K/hr calculation
  sqbot = stn.([qbotfld,'_err']).data + (0.03.*abs(qbo));	s2qbot = signedSQR(sqbot);
  stn.([qbotfld,'_err']).data(:) = sqbot(:);


  COVgQ0Qb = cov_ts(sgQ0,sQb,@isfinite,@isfinite);

  s2bQ0 = s2gQ0 + s2Qb + 2.*COVgQ0Qb;				sbQ0 = signedsqrt(s2bQ0);

  stn.([bq0fld '_err']).date = dts(:);
  stn.([bq0fld '_err']).data = sbQ0(:);


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% STORED HEAT (Q0+Qb)/rho/Cp/h

  stn = station_heat_flux_term(stn,[bq0fld,'_err'],[bq0tfld,'_err'],sfld,sal,mhfld);
  % Adjust for uncertainty of T (0.5K), S (1.0PSU), and h (18m) in K/hr calculation
  sbQ0t = stn.([bq0tfld,'_err']).data + (sKperHr.*abs(bq0)) + (5e-5.*abs(qbo));
  s2bQ0t = signedSQR(sbQ0t);
  stn.([bq0tfld,'_err']).data(:) = sbQ0t(:);


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% KM-SCALE HEAT ADVECTION

  sfqudT = 0.03.*abs(fqudT);					s2fqudT = signedSQR(sfqudT);
  stn.([fqudTfld '_err']).date = dts(:);
  stn.([fqudTfld '_err']).data = sfqudT(:);

  stn = station_heat_flux_term_inverse(stn,[fqudTffld,'_err'],[fqudTfld,'_err'],sfld,sal,mhfld);
  sfqudTf = stn.([fqudTffld,'_err']).data(:) + (sWm2.*abs(fqudT));
  stn.([fqudTffld,'_err']).data(:) = sfqudTf(:);

  COVbQ0tfqudT = cov_ts(sbQ0t,sfqudT,@isfinite,@isfinite);

  s2qtAdv = s2bQ0t + s2fqudT + 2.*COVbQ0tfqudT;			sqtAdv = signedsqrt(s2qtAdv);

  stn.([qtAdvfld '_err']).date = dts(:);
  stn.([qtAdvfld '_err']).data = sqtAdv(:);
  stn = station_heat_flux_term_inverse(stn,[qtAdvffld,'_err'],[qtAdvfld,'_err'],sfld,sal,mhfld);
  % Adjust for uncertainty of T (~0.5K), S (~1.0PSU), and h (~20m) in K/hr calculation
  sqtAdvf = stn.([qtAdvffld,'_err']).data(:) + (sWm2.*abs(fqudT)); s2qtAdvf = signedSQR(sqtAdvf);
  stn.([qtAdvffld,'_err']).data(:) = sqtAdvf(:);


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% FICKIAN HEAT DIFFUSION

  skd2T = 0.03.*abs(kd2T);					s2kd2T = signedSQR(skd2T);
  stn.([kd2Tfld '_err']).date = dts(:);
  stn.([kd2Tfld '_err']).data = skd2T(:);
  stn = station_heat_flux_term_inverse(stn,[kd2Tffld,'_err'],[kd2Tfld,'_err'],sfld,sal,mhfld);
  skd2Tf = stn.([kd2Tffld,'_err']).data(:);


  COVqtAdvkd2T = cov_ts(sqtAdv,skd2T,@isfinite,@isfinite);

  s2bdT = s2qtAdv + s2kd2T + 2.*COVqtAdvkd2T;			sbdT = signedsqrt(s2bdT);

  stn.([bdTfld '_err']).date = dts(:);
  stn.([bdTfld '_err']).data = sbdT(:);
  stn = station_heat_flux_term_inverse(stn,[bdTffld,'_err'],[bdTfld,'_err'],sfld,sal,mhfld);
  sbdTf = stn.([bdTffld,'_err']).data(:) + (sWm2.*(abs(fqudT)+abs(kd2T)));	s2bdTf = signedSQR(sbdTf);
  stn.([bdTffld,'_err']).data(:) = sbdTf(:);


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% TOTAL HEAT BUDGET

  shchc = 0.03.*abs(hchc);					s2hchc = signedSQR(shchc);

  stn.([hcdTdthc '_err']).date = dts(:);
  stn.([hcdTdthc '_err']).data = shchc(:);
  stn = station_heat_flux_term_inverse(stn,[hcdTdthcf,'_err'],[hcdTdthc,'_err'],sfld,sal,mhfld);
  shchcf = stn.([hcdTdthcf,'_err']).data(:);			s2hchcf = signedSQR(shchcf);

  COVbdThchc = cov_ts(sbdT,shchc,@isfinite,@isfinite);
  s2hc = s2bdT + s2hchc + 2.*COVbdThchc;			shc = signedsqrt(s2hc);

  stn.([hcdTdt '_err']).date = dts(:);
  stn.([hcdTdt '_err']).data = shc(:);

  %%%%DEBUG???
  %stn = station_heat_flux_term(stn,[sq0fld,'_err'],[hcdTdt,'_err'],sfld,sal,mhfld);

  stn = station_heat_flux_term_inverse(stn,[hcdTdtf,'_err'],[hcdTdt,'_err'],sfld,sal,mhfld);
  shcf = stn.([hcdTdtf,'_err']).data(:);			s2hcf = signedSQR(shcf);


  %toc,
  % timenow,
  set_more;

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% INTERNAL FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = signedSQR(val)
  val = sign(val).*(val.^2);
return;

function val = signedsqrt(val)
  val = sign(val).*sqrt(abs(val));
return;
