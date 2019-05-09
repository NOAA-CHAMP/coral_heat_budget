1;

fld = [];

if ( ~exist('bbox','var') )
  bbox = [24,27,-80,-78];
end;

if ( ~exist('begdt','var') )
  begdt = datenum(1993,1,1,0,0,0);
end;

%DEBUG:enddt = datenum(1993,1,1,18,0,0);
%DEBUG:enddt = datenum(1993,1,31,18,0,0);
if ( ~exist('enddt','var') )
  %enddt = datenum(2014,7,1,0,0,0);
  enddt = datenum(2014,6,30,18,0,0);
end;

datapath = get_thesis_path('../data');
matfpat = sprintf('erai_%06.3f_%06.3f_%06.3f_%06.3f_%s_%s',bbox,datestr(begdt),datestr(enddt));
matfpat = strrep(matfpat,':','_');
matfname = fullfile(datapath,[matfpat,'.mat']);

if ( exist(matfname,'file') )
  disp(['LOAD ',matfname]);
  load(matfname,'fld');

else
  eraipath = get_thesis_path('../data/ECMWF/ERAI');
  fnamepat = ['ei.oper.an.sfc.regn128sc.%04d%02d%02d%02d.gramer92187.nc'];

  fld.bbox = bbox;
  fld.begdt = begdt;
  fld.enddt = enddt;

  fld.nDts = (enddt-begdt)*4;
  %fld.nLats = min(32,ceil((bbox(2)-bbox(1)+2)/0.7018)); %max 32
  %fld.nLons = min(57,ceil((bbox(4)-bbox(3)+2)/0.7031)); %max 57

  dtix = 0;
  fld.dts = repmat(nan,[fld.nDts,1]);
  oldyr = nan;
  for dt=begdt:enddt %yr=1993:2014
    [yr,mo,dy,ig,ig,ig] = datevec(dt);

    %DEBUG:  disp(datestr(dt);
    %DEBUG:
    if yr ~= oldyr; oldyr = yr; disp(datestr(dt)); end;

    for hr=0:6:18
      fname = fullfile(eraipath,sprintf(fnamepat,yr,mo,dy,hr));
      if ( ~exist(fname,'file') )
        warning('MISSING FILE: %s',fname);
        continue;
      end;

      dtix = dtix + 1;
      fld.dts(dtix) = datenum(yr,mo,dy,hr,0,0);

      nc = mDataset(fname);
      if ( isempty(nc) )
        warning('INVALID FILE: %s',fname);
        clear nc;
        continue;
      end;

      try,
        if ( ~isfield(fld,'allLats') )
          fld.allLats = cast(nc{'g4_lat_1'}(:),'double');
          fld.allLons = cast(nc{'g4_lon_2'}(:),'double');
          fld.allLons(fld.allLons >= 180) = fld.allLons(fld.allLons >= 180) - 360;
          fld.latix = find(bbox(1)-1<=fld.allLats & fld.allLats<=bbox(2)+1);
          fld.lonix = find(bbox(3)-1<=fld.allLons & fld.allLons<=bbox(4)+1);

          fld.nLats = numel(fld.latix);
          fld.nLons = numel(fld.lonix);
          fld.lats = fld.allLats(fld.latix);
          fld.lons = fld.allLons(fld.lonix);

          fld.u = repmat(nan,[fld.nDts,fld.nLats,fld.nLons]);
          fld.v = repmat(nan,[fld.nDts,fld.nLats,fld.nLons]);
          fld.Ta = repmat(nan,[fld.nDts,fld.nLats,fld.nLons]);
          fld.Td = repmat(nan,[fld.nDts,fld.nLats,fld.nLons]);
          fld.A = repmat(nan,[fld.nDts,fld.nLats,fld.nLons]);
          fld.Ts = repmat(nan,[fld.nDts,fld.nLats,fld.nLons]);
          fld.C = repmat(nan,[fld.nDts,fld.nLats,fld.nLons]);
        end;

        dat = cast(nc{'10U_GDS4_SFC'}(:,:,:),'double');
        fld.u(dtix,:,:) = dat(fld.latix,fld.lonix);
        dat = cast(nc{'10V_GDS4_SFC'}(:,:,:),'double');
        fld.v(dtix,:,:) = dat(fld.latix,fld.lonix);
        dat = cast(nc{'2T_GDS4_SFC'}(:,:,:),'double');
        fld.Ta(dtix,:,:) = dat(fld.latix,fld.lonix)-273.14;
        dat = cast(nc{'2D_GDS4_SFC'}(:,:,:),'double');
        fld.Td(dtix,:,:) = dat(fld.latix,fld.lonix)-273.14;
        dat = cast(nc{'10U_GDS4_SFC'}(:,:,:),'double');
        fld.A(dtix,:,:) = dat(fld.latix,fld.lonix);
        dat = cast(nc{'SSTK_GDS4_SFC'}(:,:,:),'double');
        fld.Ts(dtix,:,:) = dat(fld.latix,fld.lonix)-273.14;
        dat = cast(nc{'TCC_GDS4_SFC'}(:,:,:),'double');
        fld.C(dtix,:,:) = dat(fld.latix,fld.lonix);

        close(nc); clear nc

      catch,
        catchwarn(mfilename);
        close(nc); clear nc
      end;

    end; %for hr

  end; %for dt

  disp(['SAVE ',matfname]);
  save(matfname,'fld');

end; %if exist
