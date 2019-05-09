function stn = query_glb_analysis_subset(stn,mindt,maxdt,vars,baseurl)
%function stn = query_glb_analysis_subset(stn,mindt,maxdt,vars,baseurl)
%
% Add fields for ocean surface U and V current components from assimilative,
% operational NRL Global HYCOM (1/12 degree) hydrodynamic model, to struct
% STN.  Expects STN to already have .lon and .lat fields with its position.
%
% Optional VARS is a string or cell array of strings specify which variables
% to extract. The catalog of available variables can be found here:
%     http://tds.hycom.org/thredds/dodsC/glb_analysis.html
%
% BASEURL defaults to http://tds.hycom.org/thredds/dodsC/glb_analysis.ascii';
% if you have another ocean model with a THREDDS interface, specify it here.
%
% CALLS: QUERY_GLB_ANALYSIS_INDICES, URLREAD, TEXTSCAN
%
% Last Saved Time-stamp: <Fri 2010-05-07 22:14:17 Eastern Daylight Time Lew.Gramer>

  zerodt = datenum(2003,11,16);

  if ( ~exist('mindt','var') || isempty(mindt) )
    mindt = zerodt;
  end;
  if ( ~exist('maxdt','var') || isempty(maxdt) )
    % Takes a day or two to post these run results sometimes
    maxdt = floor(now - 2);
  end;
  if ( ~exist('vars','var') || isempty(vars) )
    %vars = { 'u', 'v', 'temperature', 'salinity' };
    vars = { 'u', 'v' };
  end;

  if ( ~exist('baseurl','var') || isempty(baseurl) )
    baseurl = 'http://tds.hycom.org/thredds/dodsC/glb_analysis.ascii';
  end;

  % EXAMPLE: minix = 2353 for datenum(2010,04,29)
  minix = floor(mindt) - zerodt;
  maxix = floor(maxdt) - zerodt;
  ixes = minix:maxix;

  dts = mindt:maxdt;
  ndts = numel(dts);

  [yix,xix] = query_glb_analysis_indices(stn.lon,stn.lat);


  % Would love to do subset via DODS and netCDF! But none of the following work well...
  %nc = mDataset('http://tds.hycom.org/datasets/hycom/global/glb_analysis/data/2d/archv.2010_080_00_2d.nc');
  %http://tds.hycom.org/thredds/dodsC/glb_analysis.dods?Latitude[1822][2559],Longitude[1822][2559],Date[0:1:2321],mixed_layer_u_velocity[0:1:2321][1822][2559],mixed_layer_v_velocity[0:1:2321][1822][2559],mixed_layer_thickness[0:1:2321][1822][2559],mixed_layer_temperature[0:1:2321][1822][2559],mixed_layer_salinity[0:1:2321][1822][2559],Depth[0:1:8],u[0:1:2321][0:1:8][1822][2559],v[0:1:2321][0:1:8][1822][2559],temperature[0:1:2321][0:1:8][1822][2559],salinity[0:1:2321][0:1:8][1822][2559]
  %close(nc);

  % So when life offers you ASCII-only lemons... parse that ASCII!
  for vix = 1:length(vars)
    var = vars{vix};
    stn.(['global_hycom_' var]).date = dts';

    for begix = 1:300:length(dts)
      endix = min((begix + 300 - 1), length(dts));

      url = sprintf('%s?%s[%d:1:%d][0:1:0][%d][%d]', ...
                    baseurl, var, ixes(begix), ixes(endix), yix, xix);
      s = urlread(url);
      ix = strfind(s,sprintf('\n['));
      c = textscan(s(ix(1):end),'%*[[]%*d%*[^,], %f\n', ndts);
      stn.(['global_hycom_' var]).data(begix:endix,1) = c{:};
    end;
  end;


return;
