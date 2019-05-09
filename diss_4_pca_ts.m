function [pc,evecs,evls,explnd] = diss_4_pca_ts(stn,fld,mapdims,period,ndays,evecs,evls,explnd)
%function [pc,evecs,evls,explnd] = diss_4_pca_ts(stn,fld,mapdims,period,ndays,evecs,evls,explnd)

  doPrint = true;

  if ( ~exist('fld','var') || isempty(fld) )
    fld = 'ndbc_sea_t';
  end;
  if ( ~exist('mapdims','var') || isempty(mapdims) )
    mapdims = [1,3];
  end;
  if ( ~exist('period','var') || isempty(period) )
    period = 'annual';
  end;
  if ( ~exist('ndays','var') || isempty(ndays) )
    ndays = 1;
  end;

  stn = verify_variable(stn,fld);
  if ( ~isfield(stn,fld) )
    error('No such field %s!',fld);
  end;

  nvars = 1;

  nhrs = ndays*24;

  stnm = lower(stn.station_name);

  [fdts, fdat] = som_vs_pca_ts_preprocess(stn.(fld).date, stn.(fld).data,period,ndays,[],'linear');
  varmean = repmat(nanmean(fdat), [size(fdat,1) 1]); fanom = fdat - varmean;
  fanom = fdat - varmean;

  maxix = floor((365*3) / ndays);
  %maxix = floor((3653/2) / ndays);
  % maxix = floor((3653) / ndays);
  % maxix = floor((2*3653) / ndays);
  maxix = min(maxix,size(fdat,1));

  ixes = 1:maxix;
  begyr = get_year(fdts(ixes(1)));
  endyr = get_year(fdts(ixes(end)));

  clear fdts fdat varmean


  if ( ~exist('evecs','var') || numel(evecs) ~= (numel(ixes)^2) )
    warning('Rerunning PCA');
    pause(3);
    evecs=[]; evls=[]; explnd=[];
    clear evecs evls explnd
    covm = cov(fanom(ixes,:)');
    [evecs, evls, explnd] = pcacov(covm);
    covm=[]; clear covm
  end;

  pc=fanom(ixes,:)' * evecs;
  pc=pc';

  % If PCA fouls up, e.g., the diurnal warming cycle
  %pc=-pc;


  nmaps = mapdims(2)*mapdims(1);
  ncomps = min(size(pc,1), nmaps);
  totpct = 0;

  fh = fmg;

  cmins = min(pc(1:ncomps,:)); cmn = min(reshape(cmins, [nhrs nvars]));
  cmaxs = max(pc(1:ncomps,:)); cmx = max(reshape(cmaxs, [nhrs nvars]));
  cmin = min([-abs(cmn) ; -abs(cmx)]); cmax = max([abs(cmn) ; abs(cmx)]);
  cmin = cmn; cmax = cmx;
  for ix = 1:ncomps
    if ( ix > nmaps )
      break;
    end;
    subplot(mapdims(2), mapdims(1), ix);
    hold on;
    x = 1:nhrs;
    pcadats = reshape(pc(ix,:), [nhrs nvars]);
    diss_4_pca_ts_plot(x, pcadats, cmin, cmax);
    hold off;
    title(sprintf('#%d (%0.1f%% %s)', ix, explnd(ix), '\sigma^2'));
    totpct = totpct + explnd(ix);
  end;

  suplabel(strrep(lower(fld),'_','\_'), 'x');
  suplabel(sprintf( 'PCA(%s) Modes %s %d-%d %s, %d h frames (%3.1f%% of %s): %s%s', ...
                    '\rho', upper(stnm), begyr, endyr, period, ...
                    nhrs, totpct, '\sigma^2', ...
                    strrep(lower(''),'_','\_'), '' ), 't');

  if (doPrint)
    print('-dtiff',fullfile(get_thesis_path('../DISS'),[mfilename,'-',stnm,'-',fld,'-',num2str(nhrs),'h.tif']));
  end;

return;
