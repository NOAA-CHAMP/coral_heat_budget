1;
% SCRIPT ERRBDGT
%
% Calculate full heat budget hourly errors using analytical error propagation formulae
%
% Last Saved Time-stamp: <Thu 2011-01-27 08:54:09  lew.gramer>


clear Albedo               d2                   s2Ug                 
clear COVQlhQsh            datapath             s2Ug2                
clear COVQlwQlh            dix                  s2V                  
clear COVQlwQsh            doPlot               s2V2                 
clear COVU2Ug2             dts                  s2a                  
clear COVU2V2              g                    s2d                  
clear COVUUg               gamma                s2expCair            
clear COVUV                gamma2               s2expCdew            
clear COVV2Ug2             h                    s2gamma              
clear COVVUg               h0                   s2gammaQsw           
clear COVgammaQswQlh       h02                  s2gammaQsw1          
clear COVgammaQswQlw       h0ix                 s2gammaQsw2          
clear COVgammaQswQsh       h2                   s2gammaQsw3          
clear COVqsqa              hix                  s2h                  
clear COVsa                ig                   s2q                  
clear Cair                 ix                   s2qa                 
clear Cd                   q                    s2qs                 
clear Cd2                  q0                   s2s                  
clear Cdew                 q02                  s2th                 
clear Ce                   q0ix                 s2w                  
clear Ce2                  q2                   sQ0                  
clear Ch                   qa                   sQlh                 
clear Ch2                  qa2                  sQlw                 
clear Cp                   qaix                 sQlwi                
clear Cp2                  qeix                 sQlwo                
clear Cpa                  qelt                 sQsh                 
clear Cpa2                 qlh                  sQsw                 
clear Le                   qlh2                 sQswi                
clear Le2                  qlhix                sRH                  
clear ONEHOUR              qlw                  sU                   
clear RH                   qlw2                 sU2                  
clear RH2                  qlwix                sUg                  
clear RHix                 qs                   sUg2                 
clear U                    qs2                  sV                   
clear U2                   qsh                  sV2                  
clear Ug                   qsh2                 sa                   
clear Ug2                  qshix                salin                
clear Ugix                 qsix                 sd                   
clear Uix                  qsw                  sexpCair             
clear V                    qsw2                 sexpCdew             
clear V2                   qswix                sgamma               
clear Vix                  r_sRH_sexp006a       sgammaQsw            
clear a                    r_sexpCdew_sexpCair  sh                   
clear a2                   r_sgamma_sQsw        six                  
clear airdens              r_sq_sw              sq                   
clear airdens2             r_sth_sw             sqa                  
clear aix                  rho                  sqs                  
clear alph                 rho2                 srfix                
clear ans                  s                    ss                   
clear asrfix               s2                   sth                  
clear badix                s2Q0
clear bulkix               s2Qlh                sw                   
clear s2Qlw                swh                  
clear s2Qlwi               swh2                 
clear s2Qlwo               swhix                
clear s2Qsh                th                   
clear s2Qsw                th2                  
clear corix                s2Qswi               w                    
clear cors2Q0              s2RH                 w2                   
clear corsQ0               s2U                  
clear d                    s2U2

clear AC_sQ0          hiWindIx        oldres          s2uHC           suHC
clear C               hr1             qerr            s2uf            suf
clear C2              hr2             r_sH0_sh        s2uqe           suqe
clear Fc              insuf           r_suHC_sdelTHC  sH0             uf
clear Tf              corinsuf        loWindIx        r_suf_sh        sSWH
clear bet             cornewq         newq            s2H0            sdelTHC
clear bet2            cornewres       newres          s2SWH           sdiff
clear cordts          corqerr         oldq            s2delTHC        uf2


doPlot = false;

datapath = get_thesis_path('../data');


% CONSTANTS

g = 9.8;
ONEHOUR = 3600;

airdens = 1.2;
airdens2 = airdens.^2;
Cpa = 1004.7;   % heat capacity of air [J/kg/K]
Cpa2=Cpa.^2;


% CALCULATED DATA

if ( ~isfield(stn,'ndbc_dew_t') )
  if ( isfield(stn,'ndbc_relhumid') )
    stn = station_relhumid_to_dewp(stn,'ndbc_air_t','ndbc_relhumid','ndbc_dew_t');
  else
    x = load_all_ndbc_data([],'smkf1');
    stn.ndbc_dew_t = x.ndbc_dew_t;
    x=[]; clear x;
  end;
end;
if ( ~isfield(stn,'ndbc_relhumid') )
  stn = station_dewp_to_relhumid(stn,'ndbc_air_t','ndbc_dew_t','ndbc_relhumid');
end;
if ( ~isfield(stn,'ndbc_sea_spechumid') )
  stn.ndbc_sea_spechumid.date = stn.ndbc_sea_t.date;
  % 0.98 factor from Stommel - accounts for salinity
  stn.ndbc_sea_spechumid.data = 0.98 .* relhumid_to_spechumid(stn.ndbc_sea_t.data,100);
end;
if ( ~isfield(stn,'ndbc_ncep_30a_cordiags') )
  x = load(fullfile(datapath,[lower(stn.station_name) '_cordiags.mat']));
  stn.ndbc_ncep_30a_cordiags = x.station.ndbc_ncep_30a_cordiags;
  x=[]; clear x;
end;

% Actual time series of bulk coefficients from COARE 3.0a algorithm are used below
% Cd = 1.1e-3;
% Ch = 1.0e-3;
% Ce = 1.3e-3;

[six,aix,dix,RHix,qaix,qsix,Uix,Vix,Ugix,hix,swhix,qswix,qlwix,qlhix,qshix,q0ix,h0ix,qeix] = ...
    intersect_all_dates([], stn.ndbc_sea_t.date,stn.ndbc_air_t.date,stn.ndbc_dew_t.date,...
                        stn.ndbc_relhumid.date,stn.ndbc_spechumid.date,stn.ndbc_sea_spechumid.date,...
                        stn.ndbc_wind1_u.date,stn.ndbc_wind1_v.date,...
                        stn.ndbc_ncep_30a_cordiags.date,...
                        stn.tmd_tide_i_depth.date,...
                        stn.ww3_sigwavehgt.date, ...
                        stn.ncep_dsrf.date,...
                        stn.ncep_dlrf.date,...
                        stn.ndbc_ncep_30a_latent_heat_flux.date,...
                        stn.ndbc_ncep_30a_sensible_heat_flux.date,...
                        stn.ndbc_ncep_30a_absorbed_heat_flux.date,...
                        stn.ndbc_ncep_30a_absorbed_heat_flux_24_hour_average.date,...
                        stn.fkeys_hycom_qelt.date);

dts = stn.ndbc_sea_t.date(six);

s = stn.ndbc_sea_t.data(six);					s2 = s.^2;
a = stn.ndbc_air_t.data(aix);					a2 = a.^2;
d = stn.ndbc_dew_t.data(dix);					d2 = d.^2;
RH = stn.ndbc_relhumid.data(RHix);				RH2 = RH.^2;
qs = stn.ndbc_sea_spechumid.data(qsix);				qs2 = qs.^2;
qa = stn.ndbc_spechumid.data(qaix);				qa2 = qa.^2;
U = stn.ndbc_wind1_u.data(Uix);					U2 = U.^2;
V = stn.ndbc_wind1_v.data(Vix);					V2 = V.^2;
Ug = real(stn.ndbc_ncep_30a_cordiags.Ug(Ugix));			Ug2 = Ug.^2;
Cd = real(stn.ndbc_ncep_30a_cordiags.Cd(Ugix));			Cd2 = Cd.^2;
Ch = real(stn.ndbc_ncep_30a_cordiags.Ch(Ugix));			Ch2 = Ch.^2;
Ce = real(stn.ndbc_ncep_30a_cordiags.Ce(Ugix));			Ce2 = Ce.^2;
h = stn.tmd_tide_i_depth.data(hix);				h2 = h.^2;
swh = stn.ww3_sigwavehgt.data(swhix);				swh2 = swh.^2;
qsw = stn.ncep_dsrf.data(qswix);				qsw2 = qsw.^2;
qlw = stn.ncep_dlrf.data(qlwix);				qlw2 = qlw.^2;
qlh = real(stn.ndbc_ncep_30a_latent_heat_flux.data(qlhix));	qlh2 = qlh.^2;
qsh = real(stn.ndbc_ncep_30a_sensible_heat_flux.data(qshix));	qsh2 = qsh.^2;

q0 = real(stn.ndbc_ncep_30a_absorbed_heat_flux.data(q0ix));
q0(~isfinite(q0)) = 0;
q02 = q0.^2;

qelt = real(stn.fkeys_hycom_qelt.data(qeix));
qelt(~isfinite(qelt)) = 0;


% Seawater salinity
salin = repmat(36,size(s));
% Seawater density
rho = sw_dens0(salin,s);					rho2 = rho.^2;
% Heat capacity of seawater [J/kg/K], Cp ~ 4e3
Cp = sw_cp(salin,s,h);						Cp2=Cp.^2;
% Seawater thermal expansion coefficient
alph = sw_alpha(salin,s,h);

% H0 = One Hour * Buoyancy flux / (g*alpha) = One Hour * Q0 / (rho*Cp)
h0 = ONEHOUR.*real(stn.ndbc_ncep_30a_absorbed_heat_flux_24_hour_average.data(h0ix)) ./ (rho.*Cp);
h02 = h0.^2;


% ERROR ESTIMATES

% First, estimate sQlh and sQsh using COARE 3.0a methodology: THIS IS SLOW!
if ( ~exist('cordts','var') )
  [cordts,corhs,cordwth,corhl,cordwq] = corerrbdgt(stn);
end;

% MEASUREMENT ERRORS - see NDBC, 2009
ss = 0.08;		% Sea temperature
s2s = ss.^2;
sa = 0.09;		% Air temperature
s2a = sa.^2;
sd = 0.31;		% Dewpoint temperature
s2d = sd.^2;
sU = 0.55 + (0.03.*U);	% Wind U-component [m/s] (speed+/-0.55m/s, dir+/-9oT)
s2U = sU.^2;
sV = 0.55 + (0.03.*V);	% Wind V-component [m/s]
s2V = sV.^2;

% sQswi = 50;		% Insolation peak (summer) bias in NARR [Markovic et al 2009]
% sQswi = 7;		% Insolation bias in ERA40 [Markovic et al 2009]
% Consistent with Markovic et al (2009), 14% relative error in NARR Qsw
sQswi = 50 + (0.14.*qsw);					s2Qswi = sQswi.^2;

% sQlwi = 10;		% DLRF peak (winter) bias in NARR [Markovic et al 2009]
% sQlwi = 3;		% DLRF bias in ERA40 [Markovic et al 2009]
% Consistent with Markovic et al (2009), 3% relative error in NARR Qlw
sQlwi = 10 + (0.03.*qlw);					s2Qlwi = sQlwi.^2;

% ULRF error estimated by regressing QlwO vs. in situ blackbody radiation
sQlwo = 4;							s2Qlwo = sQlwo.^2;

% TMD Tide error based on SMKF1 tide sensor: regression slope 0.7, RMSE 0.1[m]
s2h = (0.1 + (0.3.*h));						sh = sqrt(abs(s2h));
% TRNF1 slope 1.01, RMSE 0.1[m], bias -1.1[m]; SMKF1 bias was 0.75[m]

% Latent heat of evaporation [J/kg], Le ~ 2.44e6: Assumed CERTAIN in calculations
Le = vapor(nanmean(stn.ndbc_air_t.data(aix)));			Le2 = Le.^2;

% Gustiness parameter 'Ug' - error approximated based on Ug ~ 6/airtemp
s2Ug = (36./(a.^4)).*s2a;					sUg = sqrt(abs(s2Ug));

s2U2 = 4.*(U.^2).*s2U;						sU2 = sqrt(abs(s2U2));
s2V2 = 4.*(V.^2).*s2V;						sV2 = sqrt(abs(s2V2));
s2Ug2 = 4.*(Ug.^2).*s2Ug;					sUg2 = sqrt(abs(s2Ug2));


% Estimate Relative and specific humidity errors

% % Maximal error estimates based on natural variance
% sqs = 1.1e-4;							s2qs = sqs.^2;
% sqa = 1.8e-4;							s2qa = sqa.^2;

s2qs = (1.8e-5).*(0.98.^2).*(0.06.^2).*exp(0.12.*s).*s2s;	sqs = sqrt(abs(s2qs));

% % Calculate exponents in Relative Humidity formula
% c1 = 17.271; c2 = 237.7;
% Cdew = (c1) ./ (d + c2);
% Cair = (c1) ./ (a + c2);

% Simplify with a constant approximation - these coefficients vary very little!
Cdew = 0.066; Cair = 0.066;

sexpCdew = Cdew.*exp((+Cdew).*d).*sd;				s2expCdew = sexpCdew.^2;
sexpCair = Cair.*exp((-Cair).*a).*sa;				s2expCair = sexpCair.^2;

[ig,r_sexpCdew_sexpCair] = cov_ts(sexpCdew,sexpCair);

s2RH = 1e4.*( (exp(-2.*Cair.*a).*s2expCdew) + (exp(+2.*Cdew.*d).*s2expCair) + ...
              (2.*exp((+Cdew).*d).*exp((-Cair).*a).*sexpCdew.*sexpCair.*r_sexpCdew_sexpCair) );
sRH = sqrt(abs(s2RH));

[ig,r_sRH_sexp006a] = cov_ts(sRH,(0.06.*exp(0.06.*a).*sa));

s2qa = (1.8e-9).*( (exp(0.12.*a).*s2RH) + (RH2.*exp(0.12.*a).*(0.06.^2).*s2a) ...
                   + (2.*RH.*exp(0.06.*a).*r_sRH_sexp006a) );
sqa = sqrt(abs(s2qa));



%%%% ESTIMATE ERROR COVARIANCES - several of these are (should be) zero

COVsa = cov_ts(ss,sa);
%DEBUG:COVsa = 0;

COVqsqa = cov_ts(sqs,sqa);
%DEBUG:COVqsqa = 0;

COVUV = cov_ts(sU,sV);
%DEBUG:COVUV = 0;
COVUUg = cov_ts(sU,sUg);
%DEBUG:COVUUg = 0;
COVVUg = cov_ts(sV,sUg);
%DEBUG:COVVUg = 0;

COVU2V2 = cov_ts(sU2,sV2);
%DEBUG:COVU2V2 = 0;
COVU2Ug2 = cov_ts(sU2,sUg2);
%DEBUG:COVU2Ug2 = 0;
COVV2Ug2 = cov_ts(sV2,sUg2);
%DEBUG:COVV2Ug2 = 0;

w = sqrt(U2 + V2 + Ug2);					w2 = w.^2;

s2w = (w.^(-2)) .* ( (U2.*s2U) + (V2.*s2V) + (Ug2.*s2Ug) + (0.5.*(COVU2Ug2 + COVV2Ug2 + COVU2V2)) );
sw = sqrt(abs(s2w));

%DEBUG:nansummary(s2w);


th = s - a;							th2 = th.^2;
s2th = s2s + s2a - (2.*COVsa);					sth = sqrt(abs(s2th));

[ig,r_sth_sw] = cov_ts(sth,sw);
%DEBUG:r_sth_sw = 0;

s2Qsh = airdens2.*Cpa2.*Cd.*Ch.*( (th2.*s2w) + (w2.*s2th) + (2.*w.*th.*sw.*sth.*r_sth_sw) );
sQsh = sqrt(abs(s2Qsh));

%DEBUG:
nansummary(sQsh);

if ( doPlot )

  figure; maxigraph; hold on; plot(stn.ndbc_ncep_30a_sensible_heat_flux.date,real(stn.ndbc_ncep_30a_sensible_heat_flux.data),'k-'); plot(dts,(sign(qsh).*abs(sQsh)),'r.'); datetick3; titlename('\sigmaQ_S_H');

  [ix1,ix2] = intersect_dates(stn.ndbc_ncep_30a_sensible_heat_flux.date,dts);
  [yr,mo,dy,hr,mn,sc] = datevec(stn.ndbc_ncep_30a_sensible_heat_flux.date(ix1));
  dy = datenum(yr,mo,dy,0,0,0);
  shfld.date = unique(dy);
  shfld.data = grpstats(real(stn.ndbc_ncep_30a_sensible_heat_flux.data(ix1)),dy);
  sherr.date = unique(dy);
  sherr.data = grpstats(sQsh(ix2),dy);
  figure; maxigraph; hold on;
  linh = [];
  linh(end+1)=plot(shfld.date,shfld.data,'r.-');
  plot(sherr.date,shfld.data-sherr.data,'r+');
  plot(sherr.date,shfld.data+sherr.data,'r+');

  linh(end+1)=plot(stn.daily_oaflux_sensible_heat_flux.date,stn.daily_oaflux_sensible_heat_flux.data,'k-');
  [ix1,ix2]=intersect_dates(stn.daily_oaflux_sensible_heat_flux.date,stn.daily_oaflux_sensible_flux_err.date);
  plot(stn.daily_oaflux_sensible_heat_flux.date(ix1),[stn.daily_oaflux_sensible_heat_flux.data(ix1)-stn.daily_oaflux_sensible_flux_err.data(ix2),stn.daily_oaflux_sensible_heat_flux.data(ix1)+stn.daily_oaflux_sensible_flux_err.data(ix2)],'k+');

  datetick3;
  legend(linh,'Q_S_H \pmerr','OAFlux Q_0 \pmerr', 'Location','Best');
  titlename('Daily Sensible Heat Flux [W/m^2]');
  xlim(dts([1 end]));

end;



q = qs - qa;							q2 = q.^2;
s2q = s2qs + s2qa - (2.*COVqsqa);				sq = sqrt(abs(s2q));

[ig,r_sq_sw] = cov_ts(sq,sw);
%DEBUG:r_sq_sw = 0;

s2Qlh = airdens2.*Le2.*Cd.*Ce.*( (q2.*s2w) + (w2.*s2q) + (2.*w.*q.*sw.*sq.*r_sq_sw) );
sQlh = sqrt(abs(s2Qlh));

%DEBUG:
nansummary(sQlh);

if ( doPlot )

  figure; maxigraph; hold on; plot(stn.ndbc_ncep_30a_latent_heat_flux.date,real(stn.ndbc_ncep_30a_latent_heat_flux.data),'k-'); plot(dts,(sign(qlh).*abs(sQlh)),'r.'); datetick3; titlename('\sigmaQ_L_H');

  [ix1,ix2] = intersect_dates(stn.ndbc_ncep_30a_latent_heat_flux.date,dts);
  [yr,mo,dy,hr,mn,sc] = datevec(stn.ndbc_ncep_30a_latent_heat_flux.date(ix1));
  dy = datenum(yr,mo,dy,0,0,0);
  lhfld.date = unique(dy);
  lhfld.data = grpstats(real(stn.ndbc_ncep_30a_latent_heat_flux.data(ix1)),dy);
  lherr.date = unique(dy);
  lherr.data = grpstats(sQlh(ix2),dy);
  figure; maxigraph; hold on;
  linh = [];
  linh(end+1)=plot(lhfld.date,lhfld.data,'r.-');
  plot(lherr.date,lhfld.data-lherr.data,'r+');
  plot(lherr.date,lhfld.data+lherr.data,'r+');

  linh(end+1)=plot(stn.daily_oaflux_latent_heat_flux.date,stn.daily_oaflux_latent_heat_flux.data,'k-');
  [ix1,ix2]=intersect_dates(stn.daily_oaflux_latent_heat_flux.date,stn.daily_oaflux_latent_flux_err.date);
  plot(stn.daily_oaflux_latent_heat_flux.date(ix1),[stn.daily_oaflux_latent_heat_flux.data(ix1)-stn.daily_oaflux_latent_flux_err.data(ix2),stn.daily_oaflux_latent_heat_flux.data(ix1)+stn.daily_oaflux_latent_flux_err.data(ix2)],'k+');

  datetick3;
  legend(linh,'Q_L_H \pmerr','OAFlux Q_0 \pmerr', 'Location','Best');
  titlename('Daily Latent Heat Flux [W/m^2]');
  xlim(dts([1 end]));

end;


% Estimate short-wave radiative flux errors

% Time series of sea surface reflectances according to NCEP NARR
Albedo = stn.ncep_usrf.data(qswix) ./ qsw;
badix = (Albedo>0.5 | ~isfinite(Albedo));
Albedo(badix) = mean(Albedo(~badix));

% HACK: Assuming Albedo is time-varying but with ZERO uncertainty
s2Qsw = ((1 - Albedo).^2) .* s2Qswi;
s2Qsw(qsw<0.1) = 0;
sQsw = sqrt(abs(s2Qsw));

%DEBUG:
nansummary(sQsw);


% Assuming gamma is time-varying but with constant 7% relative uncertainty
% See ESTIM_GAMMA_ERROR.m for the basis of this relative error estimate
[asrfix,srfix] = intersect_dates(stn.ncep_srf_absorbed.date,stn.ncep_srf.date(qswix));
gamma = stn.ncep_srf_absorbed.data(asrfix) ./ stn.ncep_srf.data(qswix(srfix));
gamma(~isfinite(gamma)) = mean(gamma(isfinite(gamma)));
gamma(stn.ncep_srf_absorbed.data(asrfix) < 0.1) = 0;
gamma2 = gamma.^2;

s2gamma = 0.07 .* gamma2;
sgamma = sqrt(abs(s2gamma));
%DEBUG:
nansummary(sgamma);


[ig,r_sgamma_sQsw] = cov_ts(sgamma,sQsw);

s2gammaQsw = ( (gamma2.*s2Qsw) + (qsw2.*s2gamma) + (2.*gamma.*qsw.*sgamma.*sQsw.*r_sgamma_sQsw) );
sgammaQsw = sqrt(abs(s2gammaQsw));

%DEBUG:
s2gammaQsw1 = ( (gamma2.*s2Qsw)                                                                  );
%DEBUG:
s2gammaQsw2 = (                   (qsw2.*s2gamma)                                                );
%DEBUG:
s2gammaQsw3 = (                                     (2.*gamma.*qsw.*sgamma.*sQsw.*r_sgamma_sQsw) );


%DEBUG:
nansummary(sgammaQsw);


% Estimate long-wave radiative flux errors

% Scalar estimates from literature, and regression vs. in situ blackbody radiation
s2Qlw = (s2Qlwi + s2Qlwo);					sQlw = sqrt(abs(s2Qlw));

%DEBUG:
nansummary(sQlw);


%%%% ESTIMATE ERROR COVARIANCES FOR FINAL HEAT BUDGET TERMS
sgammaQsw(~isfinite(sgammaQsw)) = 0;
sQlw(~isfinite(sQlw)) = 0;
sQsh(~isfinite(sQsh)) = 0;
sQlh(~isfinite(sQlh)) = 0;
COVgammaQswQlw = cov_ts(sgammaQsw,sQlw);
COVgammaQswQsh = cov_ts(sgammaQsw,sQsh);
COVgammaQswQlh = cov_ts(sgammaQsw,sQlh);
COVQlwQsh = cov_ts(sQlw,sQsh);
COVQlwQlh = cov_ts(sQlw,sQlh);
COVQlhQsh = cov_ts(sQlh,sQsh);

s2Q0 = ( s2gammaQsw + s2Qlw + s2Qlh + s2Qsh + ...
         (2.*(COVgammaQswQlw + COVgammaQswQlh + COVgammaQswQsh + COVQlwQlh + COVQlwQsh + COVQlhQsh)) );
s2Q0(~isfinite(s2Q0)) = 0;
sQ0 = sqrt(abs(s2Q0));


% Now calculate net Q error based on COARE 3.0a method for turbulent fluxes
[corix,ix] = intersect_dates(cordts,dts);
% Ensure COARE estimates match our dates: where missing, fill in bulk errors
cors2Q0 = s2Q0;
cors2Q0(ix) = ...
    ( s2gammaQsw(ix) + s2Qlw(ix) + (cordwq(corix)'.^2) + (cordwth(corix)'.^2) + ...
      (2.*(COVgammaQswQlw + COVgammaQswQlh + COVgammaQswQsh + COVQlwQlh + COVQlwQsh + COVQlhQsh)) );
cors2Q0(~isfinite(cors2Q0)) = 0;
corsQ0 = sqrt(abs(cors2Q0));


%DEBUG:
nansummary(sQ0);

if ( doPlot )

  figure; maxigraph; hold on; plot(stn.ndbc_ncep_30a_absorbed_heat_flux.date,real(stn.ndbc_ncep_30a_absorbed_heat_flux.data),'k-'); plot(dts,(sign(q0).*abs(sQ0)),'r.'); datetick3; titlename('\sigmaQ_0');

  [ix1,ix2] = intersect_dates(stn.ndbc_ncep_30a_absorbed_heat_flux.date,dts);
  [yr,mo,dy,hr,mn,sc] = datevec(stn.ndbc_ncep_30a_absorbed_heat_flux.date(ix1));
  dy = datenum(yr,mo,dy,0,0,0);
  q0fld.date = unique(dy);
  q0fld.data = grpstats(real(stn.ndbc_ncep_30a_absorbed_heat_flux.data(ix1)),dy);
  q0err.date = unique(dy);
  q0err.data = grpstats(sQ0(ix2),dy);
  figure; maxigraph; hold on;
  linh = [];
  linh(end+1)=plot(q0fld.date,q0fld.data,'r.-');
  plot(q0err.date,q0fld.data-q0err.data,'r+');
  plot(q0err.date,q0fld.data+q0err.data,'r+');

  linh(end+1)=plot(stn.daily_oaflux_net_heat_flux.date,stn.daily_oaflux_net_heat_flux.data,'k-');

  datetick3;
  legend(linh,'Q_0 \pmerr','OAFlux Q_0', 'Location','Best');
  titlename('Daily Net Heat Flux [W/m^2]');
  xlim(dts([1 end]));

end;


% Quasi-Eulerian current errors

% Wave height error estimate based on NRMSE vs. TOPEX/Poseidon (Ardhuin et al 2010)
s2SWH = 0.13.*swh;						sSWH = sqrt(abs(s2SWH));

% Estimated as 1/2 mean WW3 Peak Wave Period for 2004-2008 at all sites
Fc = 0.4;

s2uqe = repmat(0,size(w));
loWindIx = find(w<=14.5);
hiWindIx = find(w> 14.5);
s2uqe(loWindIx) = 2.5e-8.*((1.25 - (0.25.*((0.5./Fc).^1.3))).^2).*4.*U2(loWindIx);
s2uqe(hiWindIx) = 2.5e-8.*((1.25 - (0.25.*((0.5./Fc).^1.3))).^2).*(14.5.^2);
s2uqe = (s2uqe.*s2U) + ((0.025.^2).*s2SWH);
suqe = sqrt(abs(s2uqe));


% Kilometer-scale heat advection errors

% Kilometer-scale heat diffusion errors


% Horizontal Convection errors

% Sum of error auto-covariances in Q0 for lags 1 to 24 hours
AC_sQ0 = 0;
for hr1=1:23; for hr2=2:23; AC_sQ0 = AC_sQ0 + cov_ts(sQ0(1+hr2:end-hr1),sQ0(1+hr2+hr1:end)); end; end;
AC_sQ0 = 0;
s2H0 = ((ONEHOUR./(rho.*Cp)).^2).*( (s2Q0./24) + (2.*AC_sQ0) );	sH0 = sqrt(abs(s2H0));

% These are expected to always be negligible!
%DEBUG:
r_suf_sh = 0;
%DEBUG:
r_sH0_sh = 0;


Tf = 24*ONEHOUR;
C = 1.20;							C2 = C.^2;
bet = stn.ngdc_offshore_slope;					bet2 = bet.^2;
uf = (g .* alph .* abs(h0).*h).^(1/3);				uf2 = uf.^2;

s2uf = (1/9).*((alph.*g.*h.*h0).^(2/3)).*...
       ( (s2h./h2) + (s2H0./h02) + (2.*sh.*sH0.*r_sH0_sh./(h.*h0)) );
suf = sqrt(abs(s2uf));

s2uHC = 9.*C2.*bet2.*((uf./h).^6).*(Tf.^4).*...
        ( (s2uf./uf2) + (s2h./h2) - (2.*suf.*sh.*r_suf_sh./(uf.*h)) );
suHC = sqrt(abs(s2uHC));

s2delTHC = (h02./(rho2.*Cp2.*h2)).*...
    ( (s2H0./h02) + (s2h./h2) + (2.*sH0.*sh.*r_sH0_sh./(h0.*h)) );
sdelTHC = sqrt(abs(s2delTHC));

[ig,r_suHC_sdelTHC] = cov_ts(suHC,sdelTHC);

% s2AdvHC = ( (delTHC2.*s2uHC) + (uHC2.*s2delTHC) + (2.*uHC.*delTHC.*suHC.*sdelTHC.*r_suHC_sdelTHC) );



% SUMMARIZE RESULTS

% The "Benefit of the Doubt" Method:
% A very conservative error estimate: As much as hourly heat budget differs
% from observed sea temperature variability, adjust as needed to eliminate as
% much of the residual as our hourly error estimates above could accommodate!
clear sdiff oldq oldres qerr newq insuf newres
sdiff = diff(s);
oldq = qelt(1:end-1);
oldres = oldq - sdiff;

qerr = ONEHOUR .* abs(sQ0) ./ (rho.*Cp.*h);
qerr = qerr(1:end-1);
newq = sdiff;
insuf = (abs(oldres) > qerr);
newq(insuf) = oldq(insuf) - (sign(oldres(insuf)).*qerr(insuf));
newres = newq - sdiff;

corqerr = ONEHOUR .* abs(corsQ0) ./ (rho.*Cp.*h);
corqerr = corqerr(1:end-1);
cornewq = sdiff;
corinsuf = (abs(oldres) > corqerr);
cornewq(corinsuf) = oldq(corinsuf) - (sign(oldres(corinsuf)).*corqerr(corinsuf));
cornewres = cornewq - sdiff;

figure; maxigraph; hold on;
plot(dts,s,'k.-');
plot(dts(2:end),s(1:end-1)+newq,'r.');
datetick3;
legend('Sea temperature','Gramer&Mariano + best bulk error', 'Location','Best');


figure; maxigraph; hold on;
plot(dts(2:end),qerr,'b-');
plot(dts(2:end),corqerr,'r.');
datetick3;
legend('Bulk \sigma[\partial_tT]','COR3 \sigma[\partial_tT]');

figure; maxigraph; hold on;
plot(dts,s,'k.-');
plot(dts(2:end),s(1)+cumsum(oldq),'bo');
plot(dts(2:end),s(1)+cumsum(newq),'r.');
plot(dts(2:end),s(1)+cumsum(cornewq),'ys');
datetick3;
legend('Sea temperature',...
       'T_0+\Sigma(Gramer&Mariano \partial_tT)',...
       'T_0+\Sigma(Gramer&Mariano \partial_tT + best bulk err)',...
       'T_0+\Sigma(Gramer&Mariano \partial_tT + best COR3 err)',...
       'Location','Best');
