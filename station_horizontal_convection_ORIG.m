function stn = station_horizontal_convection(stn,tf,sf,hf,q0f,qvf,qf,R,cfac,wfac,ntfld,hcfld,lagoff)
%function stn = station_horizontal_convection(stn,tf,sf,hf,q0f,qvf,qf,R,cfac,wfac,ntfld,hcfld,lagoff)
%
% From physical data (field names in struct STN; salinity and/or depth may be
% constant scalars, or empty to use DEFAULTS, 36 psu, and site depth, resp.): 
%
%   tf = sea temperature [oC]
%   sf = salinity [psu]
%   hf = pressure [db] or site depth [m]
%   q0f = net surface heat flux [W m^-2]
%   R = mixing(Ri) / entrainment rate (DEFAULT: 0.11)
%
% ... call THERMAL_EXCHANGE to calculate volumetric discharge [m2/s] and
% thermal exchange [K/s] rates, Qv and Qf, and add to struct STN as fields
% STN.(QVF) and STN.(QF), resp. If any of the field names TF,SF,HF,Q0F does
% not exist in STN, attempt to create it by calling VERIFY_VARIABLE (qv.).
%
% If specified, value of field STN.(NTFLD), net heating term (Q0/rho*Cp*h) is
% also added to horizontal convection and stored in new field STN.(HCFLD). If
% LAGOFF is given and non-zero, that lag (in hours) is built into addition.
%
% CALLS: INTERSECT_ALL_DATES, THERMAL_EXCHANGE, VERIFY_VARIABLE
%
% SEE: Farrow and Patterson 1993, Sturman et al. 1999, Monismith et al. 2006,
%  Hughes and Griffiths 2008, Mao et al. 2010a and 2010b, Chubarenko 2010
%
% Last Saved Time-stamp: <Tue 2010-12-28 08:57:19  lew.gramer>

  if ( ~exist('b','var') || isempty(b) )
    stn = station_ngdc_offshore_slope(stn);
    b = stn.ngdc_offshore_slope;
  end;

  if ( ~exist('R','var') || isempty(R) )
    % Prastowo et al. 2009: Mixing efficiency with hydraulic controls ~11%
    R = 0.11;
  end;

  if ( ~exist('ntfld','var') )
    ntfld = [];
  end;
  if ( ~exist('hcfld','var') )
    hcfld = 'netqf';
  end;
  if ( ~exist('lagoff','var') )
    lagoff = 0;
  end;

  stn = verify_variable(stn,tf);
  stn = verify_variable(stn,q0f);

  if ( isempty(sf) )
    sf = 36;
  elseif ( ischar(sf) )
    stn = verify_variable(stn,sf);
  elseif ( ~isnumeric(sf) || ~isscalar(sf) )
    error('Salinity arg SF must be a fieldname, constant scalar, or empty!');
  end;

  if ( isempty(hf) )
    if ( isfield(stn,'depth') )
      hf = stn.depth;
    else
      error('Depth arg HF is empty but field STN.depth was not found!');
    end;
  elseif ( ischar(hf) )
    stn = verify_variable(stn,hf);
  elseif ( ~isnumeric(hf) || ~isscalar(hf) )
    error('Depth arg HF must be a fieldname, constant scalar, or empty!');
  end;

  % Allow salinity and depth to be constant scalars
  if ( ischar(sf) && ischar(hf) )
    [tfix,sfix,hfix,q0fix] = intersect_all_dates( [], ...
        stn.(tf).date, stn.(sf).date, stn.(hf).date, stn.(q0f).date );
    s = stn.(sf).data(sfix);
    h = stn.(hf).data(hfix);
  elseif ( ischar(sf) )
    [tfix,sfix,q0fix] = intersect_all_dates( [], ...
        stn.(tf).date, stn.(sf).date, stn.(q0f).date );
    s = stn.(sf).data(sfix);
    h = repmat(hf,size(s));
  elseif ( ischar(hf) )
    [tfix,hfix,q0fix] = intersect_all_dates( [], ...
        stn.(tf).date, stn.(hf).date, stn.(q0f).date );
    h = stn.(hf).data(hfix);
    s = repmat(sf,size(h));
  else
    [tfix,q0fix] = intersect_dates(stn.(tf).date, stn.(q0f).date);
    s = repmat(sf,size(stn.(tf).date(tfix)));
    h = repmat(hf,size(s));
  end;

  dts = stn.(tf).date(tfix);
  t = stn.(tf).data(tfix);
  q0 = stn.(q0f).data(q0fix);

  [T,Dt,B,Qv,Q] = thermal_exchange(t,s,h,q0,b,R,cfac,wfac);

  stn.(qvf).date = dts;
  stn.(qvf).data = Qv;
  stn.(qf).date = dts;
  stn.(qf).data = Q;

  qtf = [qf '_onset_time'];
  stn.(qtf).date = dts;
  stn.(qtf).data = T;

  if ( isfield(stn,ntfld) )
    [q0ix,hcix] = intersect_dates(stn.(ntfld).date, stn.(qf).date);
    if (isfield(stn,hcfld)); stn = rmfield(stn,hcfld); end;
    stn.(hcfld).date = stn.(ntfld).date(q0ix(1:end-lagoff));
    stn.(hcfld).data = stn.(ntfld).data(q0ix(1:end-lagoff)) + stn.(qf).data(hcix(1+lagoff:end));
  end;

return;
