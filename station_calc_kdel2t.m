function stn = station_calc_kdel2t(stn,K_theta,fldnm,rawkl,kl,dtfld,ltfld,interpMethod,doDebug)
%function stn = station_calc_kdel2t(stn,K_theta,fldnm,rawkl,kl,dtfld,ltfld,interpMethod,doDebug)
%
% Calculate dT/dt as the sum of net ocean surface heating (Q0/rho*h*Cpe),
% advection, and Laplacian heat diffusion (K_THETA.*DEL2(T)). K_THETA may
% be scalar, a time series struct, or a climatology (sine) cycle specifier:
% a vector of three values, [mink,maxk,peak_jday], or vector of four values,
% [mink,maxk,peak_jday,period_in_days]. DEFAULT K_THETA: 20.
%
% If optional DODEBUG is not specified or TRUE, diagnostic info is DISP'd.
%
% All args except struct STN, K_THETA, and DODEBUG are FIELD NAME STRINGS:
%
% Inputs: FLDNM=sea temperature field struct; optional DTFLD=sum of net sfc.
%  heating (Q0/rho*h*Cp) and heat advection (u*gradient(T)) terms; optional 
%  arg INTERPMETHOD is passed to INTERP_FIELD (v. for DEFAULT) on Laplacian.
%
% Outputs: RAWKL=low time-resolution heat diffusion term [*K/hr*] (e.g., for
%  Global or Gulf of Mexico HYCOM products, this has one value per day);
%  KL=hourly heat diffusion term; LTFLD=hourly diffused heat KL + DTFLD.
%
% Sample call:
% >> stn = get_fkeys_hycom('mlrf1'); % Get FKEYS HYCOM data at Molasses Reef
% >> stn = station_calc_kdel2t(stn,2.5,'fkeys_hycom_seatemp_field',...
%                              'native_fkeys_hycom_diffused_heat',...
%                              'fkeys_hycom_diffused_heat',...
%                              'fkeys_hycom_qedt','fkeys_hycom_qelt');
%
% Last Saved Time-stamp: <Wed 2012-01-11 14:40:10  Lew.Gramer>

  if (~exist('K_theta','var') || isempty(K_theta))
    % From HYCOM v2.2 documentation: TEMDF2 * DX
    %    For GoM HYCOM ~ 20 [m^2/s]
    %    For FKEYS HYCOM ~ 2.5 [m^2/s]
    K_theta = 20.0e-3 * 1e3;
  end;
  if (~exist('fldnm','var') || isempty(fldnm))
    fldnm='gom_hycom_seatemp_field';
  end;
  if (~exist('rawkl','var') || isempty(rawkl))
    rawkl='raw_gom_hycom_diffused_heat';
  end;
  if (~exist('kl','var') || isempty(kl))
    kl='gom_hycom_diffused_heat';
  end;
  if (~exist('dtfld','var') || isempty(dtfld))
    dtfld=[];
  end;
  if (~exist('ltfld','var') || isempty(ltfld))
    ltfld=[];
  end;
  if (~exist('interpMethod','var') || isempty(interpMethod))
    interpMethod='linear';
  end;

  if ( ~exist('doDebug','var') || isempty(doDebug) )
    doDebug = true;
  end;


  % Make sure we have a Laplacian
  if ( ~isfield(stn.(fldnm),'laplacian') )
    stn = calc_field_terms(stn,fldnm);
  end;

  k = build_clim_opt(K_theta,'K_theta',stn.(fldnm).date,doDebug);

  l = interp_field(stn.(fldnm).lat,stn.(fldnm).lon,stn.(fldnm).laplacian,stn.lat,stn.lon,interpMethod);
  stn.([rawkl,'_calc_laplacian']).date = stn.(fldnm).date;
  stn.([rawkl,'_calc_laplacian']).data = l;
  stn.(rawkl).date = stn.(fldnm).date;
  stn.(rawkl).data = k.*3600.*l;

  %%%% ??? DEBUG
  % kllp = [rawkl '_90_day_lowpass'];
  % if (isfield(stn,kllp)); stn = rmfield(stn,kllp); end;
  % stn = verify_variable(stn,kllp);

  stn.(kl) = interp_ts(stn.(rawkl));
  %%%% ??? DEBUG
  % stn.(kl) = stn.(kllp);

  if ( ischar(ltfld) && ischar(dtfld) && isfield(stn,dtfld) )
    stn.(ltfld) = ts_op(stn.(dtfld),stn.(kl),'+');
  elseif ( ~isempty(ltfld) || ~isempty(dtfld) )
    warning('Ecoforecasts:StationKDel2T:BadFieldNames',...
            'Invalid fieldname(s) LTFLD or DTFLD for heat budget terms!');
  end;

return;
