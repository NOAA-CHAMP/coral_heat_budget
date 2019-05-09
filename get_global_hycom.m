function stn = get_global_hycom(stn_or_stnm,mindt,maxdt,vars,flds,baseurl)
%function stn = get_global_hycom(stn_or_stnm,mindt,maxdt,vars,flds,baseurl)
%
% Add fields for ocean surface U and V current components, sea temperature,
% salinity, and cross-shore temp. gradient, from the assimilative NRL Global
% HYCOM (1/12 degree) hydrodynamic model, to struct STN.  If a station name
% string STNM (five characters, e.g., 'mlrf1') is given instead of a struct,
% or if struct has no valid .lon and .lat fields (e.g., -80.38, 25.01), calls
% GET_STATION_COORDS with station name, to retrieve station location.
%
% Optional VARS is a string or cellstr specifying which variables to extract.
% The catalog of available variables can be found here:
%     http://tds.hycom.org/thredds/dodsC/glb_analysis.html
% If VARS is specified, optional FLDS may also be specified for the corresp.
% list of struct variable names to be added to STN, e.g., 'seatemp'.
%
% BASEURL defaults to http://tds.hycom.org/thredds/dodsC/glb_analysis.ascii';
% if you have another ocean model with a THREDDS interface, specify it here.
%
% SAMPLE URL:
%  http://tds.hycom.org/thredds/dodsC/glb_analysis.ascii?Latitude[1827][2568],Longitude[1827][2568],Date[0:1:1],mixed_layer_thickness[0:1:0][1827][2568],Depth[0:1:2],temperature[0:1:0][0:1:2][1827][2568],salinity[0:1:0][0:1:2][1827][2568]
%
% CALLS: QUERY_GLB_ANALYSIS_INDICES,GET_STATION_COORDS,URLREAD,TEXTSCAN
%
% Last Saved Time-stamp: <Wed 2010-09-01 15:28:01  Lew.Gramer>

  set_more off;

  if ( ischar(stn_or_stnm) )
    stn.station_name = stn_or_stnm;
  elseif ( isstuct(stn_or_stnm) )
    stn = stn_or_stnm;
  else
    error('First arg must either be station STRUCT or station name string!');
  end;

  if ( ~isfield(stn,'lon') )
    if ( isfield(stn,'station_name') )
      [stn.lon,stn.lat,stn.depth] = get_station_coords(stn.station_name);
    else
      error('Station specified without coordinates or name!');
    end;
  end;

  if ( ~exist('vars','var') || isempty(vars) )
    % vars = { 'u', 'v', 'mixed_layer_thickness', 'temperature', 'salinity' };
    % flds = { 'u', 'v', 'mld',                   'seatemp',     'salinity' };
    vars = { 'u', 'v', 'temperature', 'salinity', 'temperature',     'temperature',     };
    flds = { 'u', 'v', 'seatemp',     'salinity', 'inshore_seatemp', 'offshore_seatemp' };
  end;
  if ( ~exist('flds','var') || isempty(flds) )
    flds = lower(vars);
  end;

  if ( ~exist('baseurl','var') || isempty(baseurl) )
    baseurl = 'http://tds.hycom.org/thredds/dodsC/glb_analysis.ascii';
  end;

  % Download ALL available dates for dataset (e.g., '2.0031103E7, ...')
  % UGLY! But with gappy data and a shifting start point, unavoidable??
  %DEBUG:
  disp('Date[:]');
  url = sprintf('%s?Date', baseurl);
  s = urlread(url);
  ix = regexp(s,'\n[12][.][9012][0-9]') + 1;
  c = textscan(s(ix(1):end),'%f,');
  dts = datenum(strcat(num2str(c{:})),'yyyymmdd');

  zerodt = dts(1);
  lastdt = dts(end);

  if ( ~exist('mindt','var') || isempty(mindt) )
    mindt = zerodt;
  end;
  if ( ~exist('maxdt','var') || isempty(maxdt) )
    maxdt = lastdt;
  end;

  ixes = find(mindt <= dts & dts <= maxdt);
  dts = dts(ixes);
  % FORTRAN (and so ugly THREDDS Web interface) is 0-centered
  ixes = ixes - 1;

  ndts = numel(dts);


  [yix,xix] = query_glb_analysis_indices(stn.lon,stn.lat);

  switch ( lower(stn.station_name) )
   case {'fwyf1','lkwf1'},
    inshore_xincrement = -1;    offshore_xincrement = +1;
    inshore_yincrement =  0;    offshore_yincrement =  0;
   case {'cryf1','mlrf1'},
    inshore_xincrement = -1;    offshore_xincrement = +1;
    inshore_yincrement = +1;    offshore_yincrement = -1;
   case {'tnrf1','smkf1','sanf1','plsf1','dryf1'},
    inshore_xincrement =  0;    offshore_xincrement =  0;
    inshore_yincrement = +1;    offshore_yincrement = -1;
   otherwise,
    error('Unrecognized station name "%s"',stn.station_name);
  end;

  % I'd love to subset via DODS and netCDF! But none of the following work:
  %nc = mDataset('http://tds.hycom.org/datasets/hycom/global/glb_analysis/data/2d/archv.2010_080_00_2d.nc');
  %http://tds.hycom.org/thredds/dodsC/glb_analysis.dods?Latitude[1822][2559],Longitude[1822][2559],Date[0:1:2321],mixed_layer_u_velocity[0:1:2321][1822][2559],mixed_layer_v_velocity[0:1:2321][1822][2559],mixed_layer_thickness[0:1:2321][1822][2559],mixed_layer_temperature[0:1:2321][1822][2559],mixed_layer_salinity[0:1:2321][1822][2559],Depth[0:1:8],u[0:1:2321][0:1:8][1822][2559],v[0:1:2321][0:1:8][1822][2559],temperature[0:1:2321][0:1:8][1822][2559],salinity[0:1:2321][0:1:8][1822][2559]
  %close(nc);

  % So when life offers you ASCII-only lemons... parse that ASCII!
  for vix = 1:length(vars)
    var = vars{vix};
    fld = flds{vix};
    %DEBUG:
    disp(fld);
    stn.(['global_hycom_' fld]).date = dts';

    myyix = yix;
    myxix = xix;

    if ( strcmp(fld,'inshore_seatemp') )
      myxix = myxix + inshore_xincrement;
      myyix = myyix + inshore_yincrement;
    end;
    if ( strcmp(fld,'offshore_seatemp') )
      myxix = myxix + offshore_xincrement;
      myyix = myyix + offshore_yincrement;
    end;

    for begix = 1:300:length(dts)
      endix = min((begix + 300 - 1), length(dts));
      %DEBUG:
      disp(begix);

      url = sprintf('%s?%s[%d:1:%d][0:1:0][%d][%d]', ...
                    baseurl, var, ixes(begix), ixes(endix), myyix, myxix);
      s = urlread(url);
      ix = strfind(s,sprintf('\n['));
      c = textscan(s(ix(1):end),'%*[[]%*d%*[^,], %f\n', ndts);
      stn.(['global_hycom_' fld]).data(begix:endix,1) = c{:};
    end;
  end;

  set_more;

return;
