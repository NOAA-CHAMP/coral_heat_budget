function [stn,gxfld,gyfld] = station_calc_udotdelt(stn,unm,vnm,fldnm,tnm,rawud,ud,htfld,dtfld,interpMethod)
%function [stn,gxfld,gyfld] = station_calc_udotdelt(stn,unm,vnm,fldnm,tnm,rawud,ud,htfld,dtfld,[interpMethod])
%
% Calculate dT/dt as the sum of net ocean surface heating (Q0/rho*h*Cpe) and
% advective heat flux terms, per eqn. (1) of Gramer and Mariano, 2011.
%
% All arguments after STN struct are FIELDNAME strings:
% Inputs: UNM=ocean current U-component time series; VNM=current V; FLDNM=sea
%  temperature field struct; TNM=time series fieldname (used to build _x,_y_
%  gradient time series fieldnames); HTFLD=net surface forcing term (in the
%  case of heating, this is Q0/rho*h*Cp); optional arg INTERPMETHOD is passed
%  to CALC_UDOTDELT (v. for DEFAULT).
% Outputs: RAWUD=low time-resolution advected heat term (e.g., for Global or
%  Gulf of Mexico HYCOM currents, this time series has one value per day);
%  UD=hourly advected heat term; DTFLD=hourly advected heat UD + HTFLD.
%  STN.(GXFLD), STN.(GYFLD) will have x- and y-gradient time series, resp.
%
% If either HTFLD or DTFLD is not given, just calculate heat advection term.
%
% Last Saved Time-stamp: <Thu 2011-08-04 16:21:07  Lew.Gramer>

  if ( ~exist('htfld','var') || isempty(htfld) )
    htfld = [];
  end;
  if ( ~exist('dtfld','var') || isempty(dtfld) )
    dtfld = [];
  end;
  if ( ~exist('interpMethod','var') || isempty(interpMethod) )
    interpMethod = [];
  end;

  [uix,vix,fldix] = intersect_all_dates([],stn.(unm).date,stn.(vnm).date,stn.(fldnm).date);

  % Calculate advected heat at native time resolution of data ([m/s]*[K/m] == [K/s])
  [udotdelT,dTdx_1,dTdy_1] = calc_udotdelt(stn.(unm).data(uix),stn.(vnm).data(vix),stn.(fldnm),fldix,...
                                           interpMethod,stn.lat,stn.lon);

  gxfld = [tnm '_x'];
  gyfld = [tnm '_y'];
  stn.(gxfld).date = stn.(fldnm).date(fldix);
  stn.(gxfld).data = dTdx_1;
  stn.(gyfld).date = stn.(fldnm).date(fldix);
  stn.(gyfld).data = dTdy_1;

  stn.(rawud).date = stn.(unm).date(uix);
  % Convert to units of [K/hr]
  stn.(rawud).data = -(3600*udotdelT)';

%%%% DEBUG???
  % stn.(rawud).data(:) = 0;


  % Spline-fit an hourly time series of advected heat to native data
  stn.(ud).date = [stn.(unm).date(1):(1/24):stn.(unm).date(end)]';
  stn.(ud).data = spline(stn.(rawud).date,stn.(rawud).data,stn.(ud).date);
  stn = filter_gaps(stn,rawud,ud);
%%%% DEBUG???
  udlp = [rawud '_90_day_lowpass'];
  if (isfield(stn,udlp)); stn = rmfield(stn,udlp); end;
  stn = verify_variable(stn,udlp);

  if ( ischar(htfld) && ischar(dtfld) && isfield(stn,htfld) )
    % [ix1,ix2] = intersect_dates(stn.(htfld).date,stn.(ud).date);
    % stn.(dtfld).date = stn.(htfld).date(ix1);
    % stn.(dtfld).data = stn.(htfld).data(ix1) + stn.(ud).data(ix2);
%%%% DEBUG???
    [ix1,ix2] = intersect_dates(stn.(htfld).date,stn.(udlp).date);
    stn.(dtfld).date = stn.(htfld).date(ix1);
    stn.(dtfld).data = stn.(htfld).data(ix1) + stn.(udlp).data(ix2);
  elseif ( ~isempty(htfld) || ~isempty(dtfld) )
    warning('Ecoforecasts:StationUDotDelT:BadFieldNames',...
            'Invalid fieldname(s) HTFLD or DTFLD for heat budget terms!');
  end;

return;
