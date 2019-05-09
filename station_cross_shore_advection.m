function [stn,xgfld,lgfld]=station_cross_shore_advection(stn,orinm,unm,vnm,fldnm,tnm,rawud,ud,htfld,dtfld,interpMethod)
%function [stn,xgfld,lgfld]=station_cross_shore_advection(stn,orinm,unm,vnm,fldnm,tnm,rawud[,ud[,htfld,dtfld[,interpMethod]]])
%
% Calculate dT/dt as the sum of net ocean surface forcing (Q0/rho*h*Cp) and
% *cross-shore* advective gradient flux terms, per eqn. (1) of Gramer, 2010.
%
% All arguments after STN struct are FIELD NAME strings:
% Inputs: ORI=local isobath orientation (usually chosen for best result in
%  cross-shore advection term); UNM=ocean current U-component time series;
%  VNM=current V; FLDNM=sea temperature field struct; TNM=time series field
%  name (used to build _x,_y_,_xshore,_lshore names); HTFLD=net forcing term
%  (e.g., for heating Q0/rho*h*Cp) - added to STN.(UD) to return net budget;
%  optional INTERPMETHOD (DEFAULT 'linear') is passed to INTERP_FIELD (v.)
%
% Outputs (fields of STN): RAWUD=low time-resolution advected field term
%  (e.g., for Gulf of Mexico HYCOM currents, time series has one value/day);
%  UD=hourly advected field term; DTFLD=hourly advected field UD + HTFLD.
%
% If either HTFLD or DTFLD is not given, just calculate advected quantity.
%
% Last Saved Time-stamp: <Thu 2011-08-04 16:21:47  Lew.Gramer>

  if ( ~exist('htfld','var') || isempty(htfld) )
    htfld = [];
  end;
  if ( ~exist('dtfld','var') || isempty(dtfld) )
    dtfld = [];
  end;
  if ( ~exist('interpMethod','var') || isempty(interpMethod) )
    interpMethod = 'linear';
  end;

  if ( ischar(orinm) )
    ori = stn.(orinm);
  elseif ( isnumeric(orinm) && isscalar(orinm) )
    ori = orinm;
  else
    error('Second arg must be numeric orientation (oT) or field name in STN!');
  end;

  % Calculate cross-shore velocity component
  [stn,xfld,lfld] = station_reorient_vectors(stn,ori,unm,vnm);
  x = stn.(xfld).data;
  l = stn.(lfld).data;


  % Calc. advected field at native time resolution ([m/s]*[K/m]->[K/s] for heat)

  if ( ~isfield(stn.(fldnm),'gradient_x') || ~isfield(stn.(fldnm),'gradient_y') )
    % Calculate x- and y- gradients and field Laplacian
    stn = calc_field_terms(stn,fldnm);
  end;

  gxfld = [tnm '_x'];
  gyfld = [tnm '_y'];

  stn.(gxfld).date = stn.(fldnm).date;
  stn.(gxfld).data = ...
      interp_field(stn.(fldnm).lat,stn.(fldnm).lon,stn.(fldnm).gradient_x,stn.lat,stn.lon,interpMethod);

  stn.(gyfld).date = stn.(fldnm).date;
  stn.(gyfld).data = ...
      interp_field(stn.(fldnm).lat,stn.(fldnm).lon,stn.(fldnm).gradient_y,stn.lat,stn.lon,interpMethod);

  % Calculate along- and cross-shore field gradients
  [stn,xgfld,lgfld] = station_reorient_vectors(stn,ori,gxfld,gyfld);
% stn = get_avhrr_weekly_field(stn,true);
% stn = station_reorient_vectors(stn,'isobath_orientation','avhrr_weekly_sst_x','avhrr_weekly_sst_y');
% xgfld = 'avhrr_weekly_sst_xshore_30_day_lowpass';
% stn = verify_variable(stn,xgfld);
  xg = stn.(xgfld).data;
  lg = stn.(lgfld).data;


  [xix,xgix] = intersect_dates(stn.(xfld).date,stn.(xgfld).date);
  stn.(rawud).date = stn.(xfld).date(xix);
  % Convert to units of [K/hr]
  stn.(rawud).data = -(3600*x(xix).*xg(xgix));

%%%% DEBUG???
  % stn.(rawud).data(:) = 0;


  % Spline-fit an hourly time series of advected quantity to native data

  if ( exist('ud','var') && ischar(ud) )
    stn.(ud).date = [stn.(unm).date(1):(1/24):stn.(unm).date(end)]';
    stn.(ud).data = spline(stn.(rawud).date,stn.(rawud).data,stn.(ud).date);
    stn = filter_gaps(stn,rawud,ud);
%%%% DEBUG???
    % udlp = ud;
    udlp = [rawud '_90_day_lowpass'];
    if (isfield(stn,udlp)); stn = rmfield(stn,udlp); end;
    stn = verify_variable(stn,udlp);

    if ( ischar(dtfld) && ischar(htfld) && isfield(stn,htfld) )
      % [ix1,ix2] = intersect_dates(stn.(htfld).date,stn.(ud).date);
      % stn.(dtfld).date = stn.(htfld).date(ix1);
      % stn.(dtfld).data = stn.(htfld).data(ix1) + stn.(ud).data(ix2);
%%%% DEBUG???
      [ix1,ix2] = intersect_dates(stn.(htfld).date,stn.(udlp).date);
      stn.(dtfld).date = stn.(htfld).date(ix1);
      stn.(dtfld).data = stn.(htfld).data(ix1) + stn.(udlp).data(ix2);
    elseif ( ~isempty(htfld) || ~isempty(dtfld) )
      warning('Ecoforecasts:StationCrossShoreAdvection:BadFieldNames',...
              'Invalid fieldname(s) HTFLD or DTFLD for heat budget terms!');
    end;
  end;

return;
