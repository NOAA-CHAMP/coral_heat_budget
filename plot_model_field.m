function plot_model_field(fld,u,v,begdt,enddt,cmn,cmx)
%function plot_model_field(fld,u,v,begdt,enddt,cmn,cmx)
%
% Do a separate CONTOURF (v.) plot for each date in 2-D model field FLD (must
% be a struct with fields .lon, .lat, .date, .field, resp.) between BEGDT and
% ENDDT (DEFAULT: *all* dates).  If time series U and V are given, also plot
% ocean current vector (v. QUIVER) at the center gridpoint on each plot using
% nearest timestamp (i.e., FLD.date, U.date, V.date *need not* be identical).
% If CMN, CMX are given, set CAXIS([CMN CMX]) (v.) for all contour plots.
%
% Last Saved Time-stamp: <Wed 2011-01-19 13:49:08  lew.gramer>

  if ( ~exist('u','var') ); u = [];  end;
  if ( ~exist('v','var') ); v = [];  end;
  if ( ~exist('begdt','var') || isempty(begdt) ); begdt = fld.date(1);   end;
  if ( ~exist('enddt','var') || isempty(enddt) ); enddt = fld.date(end); end;
  if ( ~exist('cmn','var') ); cmn = [];  end;
  if ( ~exist('cmx','var') ); cmx = [];  end;

  dtix = find(begdt <= fld.date & fld.date <= enddt);
  if ( length(dtix) > 10 )
    yn=questdlg(sprintf('This will create %d plots! Continue?',length(dtix)));
    if (~strcmpi(yn,'yes') )
      return;
    end;
  end;

  % figure; maxigraph; hold on; contourf(stn.(fldnm).lon,stn.(fldnm).lat,squeeze(stn.(fldnm).field(dtix,:,:))); colorbar; quiver(stn.lon,stn.lat,stn.gom_hycom_u.data(dtix)*(3600*24/111e3),stn.gom_hycom_v.data(dtix)*0*(3600*24/111e3),0); title(sprintf('%s %s',upper(stn.station_name),datestr(stn.(fldnm).date(dtix))));

  if ( ~isempty(u) && ~isempty(v) )
    midx = round(length(fld.lon)/2);
    midy = round(length(fld.lat)/2);
    % Conversion factor from m/s to degrees (latitude) per day
    degperday = (min(diff(fld.date(:))) * 24 * 3600)/111e3;
  end;

  for ix = dtix(:)'
    dt = fld.date(ix);
    figure;
    maxigraph;
    contourf(fld.lon,fld.lat,squeeze(fld.field(ix,:,:)));
    if ( ~isempty(cmn) && ~isempty(cmx) )
      caxis([cmn cmx]);
    end;
    colorbar;
    title(datestr(dt));

    if ( ~isempty(u) && ~isempty(v) )
      [ig,uix] = min(abs(u.date - dt));
      [ig,vix] = min(abs(v.date - dt));
      hold on;
      quiver(fld.lon(midx),fld.lat(midy),...
             u.data(uix)*degperday,v.data(vix)*degperday,...
             'k:','LineWidth',2);
      hold off;
    end;
  end;

return;
