function clm = seasonalize_ww3_region(locn_or_rgn,doSeaFigs,doMonFigs,doPrint,figspath)
%function clm = seasonalize_ww3_region(locn_or_rgn,doSeaFigs,doMonFigs,doPrint,figspath)
%
% Climatologize NOAA Wave Watch III output surrounding a region RGN (SFL,
% SEF, BSC, FKS, WPR, EPR, STTSTJ, STX, or...) If DOPRINT (DEFAULT: false)
% print a PNG of each seasonal climatology in FIGSPATH (DEFAULT: directory
% heat_budget/figs). If 1st arg a 2-vector [LON,LAT] or STRUCT with fields
% .lon, .lat, optional .station_name, select RGN automatically for that loc.
%
% Last Saved Time-stamp: <Thu 2018-08-16 17:06:03 Eastern Daylight Time gramer>

  clm = [];

  rgns = {
  % Atlantic WW3 higher-res (4 m): use tight boxes to save memory

  % Florida
      'wpb',	'fknms',	[ -80.25, -79.75, 26.45, 26.95];
      'pvg',	'fknms',	[ -80.35, -79.85, 25.85, 26.35];
      'pom',	'fknms',	[ -80.40, -79.90, 25.50, 26.00];
      'uks',	'fknms',	[ -80.60, -79.90, 24.70, 25.30];
      'mks',	'fknms',	[ -81.30, -80.50, 24.40, 25.10];
      'lks',	'fknms',	[ -82.00, -81.20, 24.30, 25.00];
      % These regions are just a little too big to show detail
      'sef',	'fknms',	[ -80.25, -79.25, 25.25, 27.75];
      'bsc',	'fknms',	[ -80.75, -79.75, 25.00, 26.00];
      'fks',	'fknms',	[ -81.75, -79.75, 24.25, 25.25];
      % This region uses a large amount of memory
      'sfl',	'fknms',	[ -83.50, -79.50, 24.00, 25.50];
      % This region is too big to fit in memory!
      'bigsfl', 'fknms',	[ -83.50, -79.50, 24.00, 27.00];

  % East Caribbean
      'epr',	'pr_vi',	[ -66.00, -65.00, 17.75, 18.75];
      'wpr',	'pr_vi',	[ -67.50, -66.50, 17.75, 18.75];
      % This region uses a large amount of memory
      'pr',	'pr_vi',	[ -67.50, -65.00, 17.50, 19.00];

      'sttstj',	'pr_vi',	[ -65.10, -64.60, 18.10, 18.60];
      'stx',	'pr_vi',	[ -65.00, -64.50, 17.50, 18.00];

  % Atlantic
      'ber',	'berm',		[ -66.00, -64.00, 31.50, 33.50];

  % Pacific Islands WW3 lower-res (10 m): use big boxes for a wider view

  % Marianas
      'sai',	'cnmi', 	[ 145.50, 146.00, 14.90, 15.40];
      'tal',	'cnmi', 	[ 144.95, 145.45, 13.85, 14.35];

  % Samoas
      'pag',	'amsam', 	[-171.20,-170.20,-14.80,-13.80];
      'ofu',	'amsam', 	[-170.15,-169.15,-14.70,-13.70];
         };

  if ( ~exist('locn_or_rgn','var') || isempty(locn_or_rgn) )
    locn_or_rgn = 'stx';
  end;
  if ( ~exist('doSeaFigs','var') || isempty(doSeaFigs) )
    doSeaFigs = true;
  end;
  if ( ~exist('doMonFigs','var') || isempty(doMonFigs) )
    doMonFigs = false;
  end;
  if ( ~exist('doPrint','var') || isempty(doPrint) )
    doPrint = false;
  end;
  if ( ~exist('figspath','var') || isempty(figspath) )
    figspath = get_heat_budget_path('../figs');
  end;

  if ( ischar(locn_or_rgn) )
    rgn = lower(locn_or_rgn);
    rgnix = find(strcmpi(rgn,rgns(:,1)));
    if ( isempty(rgnix) )
      error('Do not know region "%s"',rgn);
    end;
    clm.ww3rgn = rgns{rgnix,2};
    clm.bbox = rgns{rgnix,3};

    stn.lon = mean(clm.bbox(1:2));
    stn.lat = mean(clm.bbox(3:4));
  else
    if ( isfield(locn_or_rgn,'lon') && isfield(locn_or_rgn,'lat') )
      stn = locn_or_rgn;
    elseif ( isnumeric(locn_or_rgn) && numel(locn_or_rgn) == 2 )
      stn.lon = locn_or_rgn(1);
      stn.lat = locn_or_rgn(2);
    else
      error('First arg: region name, location 2-vector, or station STRUCT');
    end;
    rgn = [];
    for rgnix=1:size(rgns,1)
      if ( ~isempty(bboxinside(stn.lon,stn.lat,rgns{rgnix,3})) )
        rgn = rgns{rgnix,1};
        clm.ww3rgn = rgns{rgnix,2};
        clm.bbox = rgns{rgnix,3};
        %DEBUG:
        disp(['Using region ',rgn]);
        break;
      end;
    end;
    if ( isempty(rgn) )
      error('Unknown region for LON %g, LAT %g',stn.lon,stn.lat);
    end;
  end;
  if ( ~isfield(stn,'station_name') )
    stn.station_name = rgn;
  end;
  clear locn_or_rgn

  clm.rgn = rgn;

  x = get_ww3_multi_region([],clm.ww3rgn,[],clm.bbox);
  clm.lon = x.lon;
  clm.lat = x.lat;
  clm.tp = x.tp;
  clm.dp = x.dp;
  clm.hs = x.hs;
  x=[]; clear x

  clm.tp.data = interp_field(clm.lat,clm.lon,clm.tp.field,stn.lat,stn.lon,@nanmean);
  clm.dp.data = interp_field(clm.lat,clm.lon,clm.dp.field,stn.lat,stn.lon,@nanmean);
  clm.hs.data = interp_field(clm.lat,clm.lon,clm.hs.field,stn.lat,stn.lon,@nanmean);

  % Surface-wave phase speed (deep-water approximation): gT/2pi
  % (If more simplicity is desired, wave period can serve as "speed".)
  spd = (9.8.*clm.tp.field)./(2*pi);
  [clm.u,clm.v] = spddir_to_uv(spd,clm.dp.field);

  maxhs = ceiln(nanmax(clm.hs.field(:)),-1);

  if ( doSeaFigs )
    for seas=1:4;
      seasix=find(get_season(clm.tp.date)==seas);
      fmg;
      contourf(clm.lon,clm.lat,squeeze(nanmedian(clm.hs.field(seasix,:,:))),[0.00:0.10:maxhs]);
      daspect([1,cosd(clm.lat(1)),1]);
      set(gca,'CLim',[0,1.5]);
      cbh=colorbar;
      ylabel(cbh,'Wave height (H_s)');
      qh = quiver(clm.lon,clm.lat,...
                  squeeze(nanmedian(clm.u(seasix,:,:))),...
                  squeeze(nanmedian(clm.v(seasix,:,:))),0.5,'k');
      %legend('H_s','c_p');
      titlename([upper(rgn),' Q-',num2str(seas),' climatology: H_s and c_p']);
      if ( isfield(stn,'ngdc_hires_bathy') )
        plot_hires_coastline(stn.ngdc_hires_bathy);
      end;
      if ( exist('set_surf_cursor','file') )
        set_surf_cursor;
      end;
      if ( doPrint )
        print('-dpng',fullfile(figspath,sprintf('ww3-hs-%s-season-%02d.png',rgn,seas)));
      end;
    end;
  end;

  if ( doMonFigs )
    for mo=1:12;
      moix=find(get_month(clm.tp.date)==mo);
      fmg;
      contourf(clm.lon,clm.lat,squeeze(nanmedian(clm.hs.field(moix,:,:))),[0.00:0.10:maxhs]);
      daspect([1,cosd(clm.lat(1)),1]);
      set(gca,'CLim',[0,1.5]);
      cbh=colorbar;
      ylabel(cbh,'Wave height (H_s)');
      qh = quiver(clm.lon,clm.lat,...
                  squeeze(nanmedian(clm.u(moix,:,:))),...
                  squeeze(nanmedian(clm.v(moix,:,:))),0.5,'k');
      %legend('H_s','c_p');
      titlename([upper(rgn),' Q-',num2str(mo),' climatology: H_s and c_p']);
      if ( isfield(stn,'ngdc_hires_bathy') )
        plot_hires_coastline(stn.ngdc_hires_bathy);
      end;
      if ( doPrint )
        print('-dpng',fullfile(figspath,sprintf('ww3-hs-%s-month-%02d.png',rgn,mo)));
      end;
    end;
  end;

return;
