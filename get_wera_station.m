function stn = get_wera_station(stn_or_stnm,xrad,yrad,interpMethod,doDownload,doWarn)
%function stn = get_wera_station(stn_or_stnm,xrad,yrad,interpMethod,doDownload,doWarn)
%
% Call ANWERA (v.) to load (2*YRAD+1)x(2*XRAD+1) time series fields (DEFAULT:
% XRAD=6, YRAD=5, 11x13 fields) of U and V surface ocean current components.
% Returns struct STN with time series 2D fields: STN.adcp_u_field.date,
% .adcp_u_field.lon, .adcp_u_field.lat, .adcp_u_field.field, etc. Then calls
% INTERP_FIELD (v.) with INTERPMETHOD (DEFAULT: 'linear') to build the time
% series STN.adcp_u and STN.adcp_v. Saves result in MAT a file, which is then
% loaded on successive calls to the function for that same STN.station_name.
%
% If optional DODOWNLOAD is *true*, ANWERA tries to download raw WERA files
% from "http://iwave.rsmas.miami.edu/wera"; otherwise, only use local files.
%
% If optional DOWARN is *false*, do not display any WERA warning messages.
%
% Last Saved Time-stamp: <Wed 2011-09-07 15:56:35  lew.gramer>

  set_more off;

  datapath = get_thesis_path('../data');

  if ( ~exist('stn_or_stnm','var') || isempty(stn_or_stnm) )
    stn_or_stnm = 'AOAT_BROAD_KEY_2';
  end;

  stn = get_station_from_station_name(stn_or_stnm);
  disp(stn.station_name);

  if ( ~exist('xrad','var') || isempty(xrad) )
    xrad = 6;
  end;
  if ( ~exist('yrad','var') || isempty(yrad) )
    yrad = 5;
  end;
  interpMethodDefaulted = false;
  if ( ~exist('interpMethod','var') || isempty(interpMethod) )
    interpMethod = 'linear';
    interpMethodDefaulted = true;
  end;
  if ( ~exist('doDownload','var') || isempty(doDownload) )
    doDownload = false;
  end;

  if ( ~exist('doWarn','var') || isempty(doWarn) )
    doWarn = true;
  end;

  prevWarns = repmat(struct('identifier',[],'state',[]),[0 0]);
  if ( ~doWarn )
    disp('Component warnings temporarily off');
    prevWarns(end+1) = warning('off','get_wera_station:BadFileShape');
    prevWarns(end+1) = warning('off','get_wera_station:SizeMismatch');
    prevWarns(end+1) = warning('off','anwera:NoURL');
    prevWarns(end+1) = warning('off','anwera:BadCoords');
    prevWarns(end+1) = warning('off','anwera:NoLocalAver');
    prevWarns(end+1) = warning('off','anwera:NoLocalAcc');
    prevWarns(end+1) = warning('off','anwera:NoData');
    prevWarns(end+1) = warning('off','anwera:BadData');
  end;


  matfname = fullfile(datapath,[lower(stn.station_name) '_wera.mat']);


  if ( exist(matfname,'file') )

    disp(['Loading ' matfname]);
    load(matfname,'result');


  else

    disp(['Processing raw WERA data files']);

    result.station_name = stn.station_name;
    result.lon = stn.lon;
    result.lat = stn.lat;

    [LONS,LATS] = anwera('20050010000',0,[],doDownload);

    nrows = size(LONS,1);
    ncols = size(LONS,2);

    [result.wera_pos_err,ix] = min(distance_wgs84(result.lat,result.lon,LATS(:),LONS(:)));

    result.wera_ix = ix;

    % Make field slightly larger on one axis, to help avoid "dyslexia errors"
    [ri,ci] = ind2sub(size(LONS),ix);
    rix = ri-yrad:ri+yrad;
    rix(1>rix | rix>nrows) = [];
    cix = ci-xrad:ci+xrad;
    cix(1>cix | cix>ncols) = [];

    nul = struct('date',[],'data',[]);
    result.wera_u = nul;
    result.wera_v = nul;
    result.wera_speed = nul;
    result.wera_dir = nul;

    result.wera_u_field.lon = unique(LONS(rix,cix));
    result.wera_u_field.lat = unique(LATS(rix,cix));
    result.wera_u_field.date = [];
    result.wera_u_field.field = repmat(nan,[0 size(LONS(rix,cix))]);

    result.wera_v_field.lon = unique(LONS(rix,cix));
    result.wera_v_field.lat = unique(LATS(rix,cix));
    result.wera_v_field.date = [];
    result.wera_v_field.field = repmat(nan,[0 size(LONS(rix,cix))]);

    result.wera_interp_method = interpMethod;

    for yr=2005:2011

      switch (yr),
       case {2004,2008,2012,2016},
        jds = 1:366;
       case 2011,
        % Data service suspended
        jds = 1:33;
       otherwise,
        jds = 1:365;
      end;

      for jd=jds(:)'

        %DEBUG:
        if (mod(jd,90)==1); disp({yr,jd}); end;

        for hr=0:23
          ds=sprintf('%04d%03d%02d00',yr,jd,hr);
          [LONS,LATS,U,V] = anwera(ds,0,[],doDownload);

          if ( ~isempty(U) )
            % CONVERT from cm/s -> m/s
            U = U./1e2;
            V = V./1e2;

            if ( nrows~=size(LONS,1) || ncols~=size(LONS,2) )
              warning('get_wera_station:BadFileShape',...
                      'ANWERA returned wrong shape for %s',ds);
            else
              dt = datenum(yr,1,1,hr,0,0) + jd - 1;
              result.wera_u_field.date(end+1,1) = dt;
              result.wera_u_field.field(end+1,:,:) = U(rix,cix);
              result.wera_v_field.date(end+1,1) = dt;
              result.wera_v_field.field(end+1,:,:) = V(rix,cix);
            end;
          end;
          clear LONS LATS U V
        end; %for hr=0:23

      end; %for jd=jds(:)'

    end; %for yr=2005:2011

    result.wera_u.date = result.wera_u_field.date;
    result.wera_u.data = interp_field(result.wera_u_field.lat, ...
                                      result.wera_u_field.lon, ...
                                      result.wera_u_field.field,...
                                      result.lat,result.lon,interpMethod);

    result.wera_v.date = result.wera_v_field.date;
    result.wera_v.data = interp_field(result.wera_v_field.lat, ...
                                      result.wera_v_field.lon, ...
                                      result.wera_v_field.field,...
                                      result.lat,result.lon,interpMethod);

    result.wera_speed.date = result.wera_u.date;
    result.wera_speed.data = uv_to_spd(result.wera_u.data,result.wera_v.data);
    result.wera_dir.date = result.wera_u.date;
    result.wera_dir.data = uv_to_dir_curr(result.wera_u.data,result.wera_v.data);


    % Calculate data density fields
    for rix=1:size(result.wera_u_field.field,2)
      for cix=1:size(result.wera_u_field.field,3)
        result.wera_u_field.n(rix,cix) = numel(find(isfinite(result.wera_u_field.field(:,rix,cix))));
        result.wera_u_field.pct(rix,cix) = result.wera_u_field.n(rix,cix)/size(result.wera_u_field.date,1);
        result.wera_v_field.n(rix,cix) = numel(find(isfinite(result.wera_v_field.field(:,rix,cix))));
        result.wera_v_field.pct(rix,cix) = result.wera_v_field.n(rix,cix)/size(result.wera_v_field.date,1);
      end;
    end;


    disp(['Saving ' matfname]);
    save(matfname,'result');


  end; %if ( exist(matfname,'file') )


  flds = fieldnames(result);
  for fldix=1:length(flds)
    fld = flds{fldix};
    stn.(fld) = result.(fld);
  end;
  result = []; clear result;

  % Just warn user if they requested a field different from what we saved
  if ( size(stn.wera_u_field.field,2) ~= ((2*yrad)+1) || ...
       size(stn.wera_u_field.field,3) ~= ((2*xrad)+1) )
    warning('get_wera_station:SizeMismatch',...
            'Requested field %dx%d does not match saved field %dx%d!',...
            ((2*yrad)+1),((2*xrad)+1),size(stn.wera_u_field.field,2),size(stn.wera_u_field.field,3));
  end;

  % Let user specify alternate interpolation methods for time series
  if ( ~strcmp(stn.wera_interp_method,interpMethod) && ~interpMethodDefaulted )
    disp(['Reinterpolating time series with ' interpMethod]);
    stn.wera_u.data = interp_field(stn.wera_u_field.lat, ...
                                      stn.wera_u_field.lon, ...
                                      stn.wera_u_field.field,...
                                      stn.lat,stn.lon,interpMethod);
    stn.wera_v.data = interp_field(stn.wera_v_field.lat, ...
                                      stn.wera_v_field.lon, ...
                                      stn.wera_v_field.field,...
                                      stn.lat,stn.lon,interpMethod);
    stn.wera_speed.data = uv_to_spd(stn.wera_u.data,stn.wera_v.data);
    stn.wera_dir.data = uv_to_dir_curr(stn.wera_u.data,stn.wera_v.data);
  end;

  if ( ~isempty(prevWarns) )
    for ix = 1:length(prevWarns)
      warning(prevWarns(ix).state,prevWarns(ix).identifier);
    end;
  end;

  set_more;

return;
