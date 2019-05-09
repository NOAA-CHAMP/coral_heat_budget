function res = antsg
%function res = antsg
%
% Return all TSG data parsed from AOML SFP (SFERP) cruises 1995-2010. Return
% struct RES will have identical length vector fields .yrs, .mos, .lon, .lat,
% .sst, .sss, .chl  from all TSG XLS files, and .d (point-to-point distance
% [km] calculated from .lon,.lat vectors), .dsst (diff(SST)/d), .keysdsst
% (DSST with all TSG points outside the Florida reef tract NaN), .crestdsst
% (all points outside FRT reef crest NaN), .z (depth from 3-arcsecond NGDC
% Coastal Relief model), .dz (topo slope), .hdg (instantaneous ship heading).
%
% Last Saved Time-stamp: <Thu 2011-05-12 10:56:01  lew.gramer>


  set_more off;

  datapath = get_thesis_path('../data');
  tsgpath = fullfile(datapath,'sfp','tsg');

  res = [];

  begyr = 1995;
  endyr = 2010;

  matfname = fullfile(datapath,'sfp_tsg.mat');
  if ( exist(matfname,'file') )

    disp(['Reloading ' matfname]);
    load(matfname,'res');

  else

    disp('Loading raw TSG data from SFP XLS files...');

    res = load('frt_reef_crest.mat');


    maxpts = round((endyr-begyr+1)*6*2500/0.3);

    res.date=repmat(nan,[maxpts,1]);
    res.yrs=repmat(nan,[maxpts,1]);
    res.mos=repmat(nan,[maxpts,1]);
    res.lon=repmat(nan,[maxpts,1]);
    res.lat=repmat(nan,[maxpts,1]);
    res.sst=repmat(nan,[maxpts,1]);
    res.sss=repmat(nan,[maxpts,1]);
    res.chl=repmat(nan,[maxpts,1]);
    res.cruise_n=repmat(nan,[maxpts,1]);

    res.hdg=repmat(nan,[maxpts,1]);
    res.d=repmat(nan,[maxpts,1]);
    res.dsst=repmat(nan,[maxpts,1]);
    res.keysdsst=repmat(nan,[maxpts,1]);
    res.crestdsst=repmat(nan,[maxpts,1]);


    ngot = 0;
    for yr=begyr:endyr
      for mo=1:12
        fpatt = fullfile(tsgpath,sprintf('%04d_%02d_*.xls',yr,mo));
        fdir = dir(fpatt);
        if ( numel(fdir) > 1 )
          fdir = fdir(1);
          warning('Duplicate cruise files: Using "%s"!',fdir.name);
        end;
        if ( ~isempty(fdir) )
          fname = fullfile(tsgpath,fdir.name);
          disp(fname);
          x = importdata(fname);

          npts = size(x.data,1);

          %DEBUG:          disp([yr,mo,npts]); continue;
          %DEBUG:          if (isfield(x,'textdata')&&size(x.textdata,2)>=6); disp({yr,mo,x.textdata{1,:}}); end; continue; % [1 2 5 6]

          lonix = find(strncmpi(x.textdata(1,:),'lon',3),1);
          latix = find(strncmpi(x.textdata(1,:),'lat',3),1);
          sstix = find((strncmpi(x.textdata(1,:),'sst',3)|strncmpi(x.textdata(1,:),'tem',3)),1);
          sssix = find((strncmpi(x.textdata(1,:),'sss',3)|strncmpi(x.textdata(1,:),'sal',3)),1);
          chlix = find((strncmpi(x.textdata(1,:),'fl',2)|strncmpi(x.textdata(1,:),'chl',3)),1);

          if ( isempty(lonix) || isempty(latix) || isempty(sstix) )
            warning('Missing lon, lat, or SST header in "%s"',fname);
            continue;
          end;

          % Date/Year/JulianDay formats are all over the place - just save the
          % year and month from the individual filenames for seasonal analysis

          lon=x.data(:,lonix);
          lat=x.data(:,latix);
          sst=x.data(:,sstix);
          if ( isempty(sssix) )
            warning('No SSS data in "%s"',fname);
            sss = repmat(nan,size(sst));
          else
            sss=x.data(:,sssix);
          end;
          if ( isempty(chlix) )
            warning('No Chl data in "%s"',fname);
            chl = repmat(nan,size(sst));
          else
            chl=x.data(:,chlix);
          end;

          % For now, assume every cruise starts on the first of the month,
          % and TSG sampling period has always and ever been one minute.
          dts = datenum(yr,mo,1) + ([1:length(lon)]'./(24*60));

          badix = find(~isfinite(lat) | ~isfinite(lon) | ~isfinite(sst));
          if ( ~isempty(badix) )
            warning('Removing %d bad (lat/lon/sst) points',numel(badix));
            dts(badix) = [];
            lon(badix) = [];
            lat(badix) = [];
            sst(badix) = [];
            sss(badix) = [];
            chl(badix) = [];
            npts = npts - numel(badix);
          end;


          [ds,hdgs] = sw_dist(lat,lon,'km');
          % Convert squirrelly SW_DIST heading -> degrees True
          % SW_DIST outputs:
          %   phaseangle  = angle of line between stations with x axis (East).
          %                 Range of values are -180..+180. (E=0, N=90, S=-90)
          hdgs = (-hdgs + 90);
          hdgs(hdgs < 0) = hdgs(hdgs < 0) + 360;

          hdg = [ nan ; hdgs ];
          d = [ nan ; ds ];

          % First value will be NaN
          dsst = [ nan ; diff(sst) ] ./ d;

          keysdsst = dsst;
          keysdsst(-81.9>lon|lon>-80.0|24.2>lat|lat>26.0) = nan;

          indix=inside(lon,lat,res.crestlon,res.crestlat)';
          crestdsst = keysdsst;
          crestdsst(indix==0)=nan;


          res.date(ngot+1:ngot+npts)=dts;
          res.yrs(ngot+1:ngot+npts)=yr;
          res.mos(ngot+1:ngot+npts)=mo;
          res.lon(ngot+1:ngot+npts)=lon;
          res.lat(ngot+1:ngot+npts)=lat;
          res.sst(ngot+1:ngot+npts)=sst;
          res.sss(ngot+1:ngot+npts)=sss;
          res.chl(ngot+1:ngot+npts)=chl;
          res.cruise_n(ngot+1:ngot+npts)=1:npts;

          res.hdg(ngot+1:ngot+npts) = hdg(:);
          res.d(ngot+1:ngot+npts) = d(:);
          res.dsst(ngot+1:ngot+npts) = dsst(:);
          res.keysdsst(ngot+1:ngot+npts) = keysdsst(:);
          res.crestdsst(ngot+1:ngot+npts) = crestdsst(:);

          ngot = ngot + npts;

          x=[]; clear x;
          clear dts lat lon sst sss chl d dsst keysdsst crestdsst;

        end; %if ( ~isempty(fdir) )
      end; %for mo=1:12
    end; %for yr=begyr:endyr

    % Trim extraneous trailing elements - if any
    flds = fieldnames(res);
    for fldix=1:length(flds)
      fld=flds{fldix};
      if ( numel(res.(fld)) >= ngot )
        res.(fld) = res.(fld)(1:ngot);
      end;
    end;


    % Not enough memory to load all 92m-resolution FRT topography at once!
    % Instead, stitch together from 81x81 km squares around each FRT site.
    stnms={'lkwf1','fwyf1','cryf1','mlrf1','lonf1','tnrf1','smkf1','looe1','amsf1','sanf1','plsf1'};
    for ix=1:length(stnms)
      stn = get_ngdc_bathy_station(stnms{ix});
      az{ix} = interp_field(unique(stn.ngdc_92m_bathy.lon),unique(stn.ngdc_92m_bathy.lat),...
                            stn.ngdc_92m_bathy.field,res.lon,res.lat,'nearest');
      stn=[]; clear stn;
    end;
    res.z = repmat(nan,size(res.lon));
    for ix=1:length(az);
      res.z(isfinite(az{ix})) = az{ix}(isfinite(az{ix}));
    end;
    az = []; clear az;

    res.dz = [ 0 ; diff(res.z) ] ./ (res.d*1e3);


    disp(['Saving MAT file ' matfname]);
    save(matfname,'res');

  end; %if ( exist(matfname,'file') ) else

  set_more;

return;


%% VERY SLOW!
% bthy=load(fullfile(get_ecoforecasts_path('coast'),'LGramer1-80.mat'));
% tic; az=repmat(nan,size(lon)); for ix=1:length(lon); [ig,zix]=min(((lon(ix)-bthy.lon).^2)+((lat(ix)-bthy.lat).^2)); az(ix)=bthy.depth(zix); end; toc,
% save('SFP-2008-Aug-transect-NGDC-depths.mat','az');
%% But already done...
az=[]; clear az;
load('SFP-2008-Aug-transect-NGDC-depths.mat','az');

keysdsst=dsst;
keysdsst(-81.9>lon|lon>-80.0|24.2>lat|lat>26.0|jd>220.5) = nan;

% indix=inside(lon,lat,cs(1,coordix),cs(2,coordix))';
indix=inside(lon,lat,crestlon,crestlat)';

crestdsst = dsst;
crestdsst(indix==0|jd>220.5)=nan;

mind=0.05;

nanmean(abs(keysdsst(mind<d&d<0.5))), nanstd(abs(keysdsst(mind<d&d<0.5))), numel(find(isfinite(keysdsst(mind<d&d<0.5)))), prctile(keysdsst(mind<d&d<0.5),[5,95]),

nanmean(abs(crestdsst(mind<d&d<0.5))), nanstd(abs(crestdsst(mind<d&d<0.5))), numel(find(isfinite(crestdsst(mind<d&d<0.5)))), prctile(crestdsst(mind<d&d<0.5),[5,95]),


mind=0.02;

nanmean(abs(keysdsst(mind<d&d<0.5))), nanstd(abs(keysdsst(mind<d&d<0.5))), numel(find(isfinite(keysdsst(mind<d&d<0.5)))), prctile(keysdsst(mind<d&d<0.5),[5,95]),

nanmean(abs(crestdsst(mind<d&d<0.5))), nanstd(abs(crestdsst(mind<d&d<0.5))), numel(find(isfinite(crestdsst(mind<d&d<0.5)))), prctile(crestdsst(mind<d&d<0.5),[5,95]),


minz=-80;

mysst=sst; mysst(jd>220.5|az<minz)=nan;
myaz=az; myaz(jd>220.5|az<minz)=nan;

fmg; [ax,h1,h2]=plotyy(cumsum(d(200:3000)),mysst(201:3001),cumsum(d(200:3000)),myaz(201:3001)); legend([h1,h2],'TSG T [^oC]','NGDC h [m]'); xlabel('Cum track dist [km]'); titlename(['SFP TSG ' datestr(datenum(yr(1),1,0)+floor(jd(1))) '-' datestr(datenum(yr(end),1,0)+floor(jd(end)))]); ylim([30,32.5]);

mydsst = diff(mysst)./d;
% fmg; boxplot(keysdsst,roundn(az(2:end),1),'notch','on');
mydaz=abs(diff(myaz))./(d*1e3);
fmg; boxplot(abs(mydsst(abs(mydsst)<2&mydaz<0.08)),roundn(mydaz(abs(mydsst)<2&mydaz<0.08),-2),'notch','on');

fmg; boxplot(abs(mydsst(abs(mydsst)<1.7&mydaz<0.1)),min(round(mydaz(abs(mydsst)<1.7&mydaz<0.1)./0.01).*0.01,0.05),'notch','on','labels',{'0','0.01','0.02','0.03','0.04','>=0.05'}); xlabel('\beta=\Deltaz/\Deltax'); ylabel('\DeltaT/\Deltax');
