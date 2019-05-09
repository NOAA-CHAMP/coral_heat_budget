function clims = anclim(stn,fld,fldix,robustAccumulator,accumulator,doPlot)
%function clims = anclim(stn,fld,fldix,robustAccumulator,accumulator,doPlot)
%
% Calculate hourly, daily, weekly, and monthly climatologies for STN.(FLD),
% and anomalies against each of the climatologies. Returned struct CLIMS has
% fields .YRS, .DTSVEC, .DATMTX, .HRCLIM, .DYCLIM, .WKCLIM, .MOCLIM, .DYMEAN,
% .WKMEAN, .MOMEAN, .HRANOM, .DYANOM, .WKANOM, and .MOANOM. FLDIX: optional
% index vector for STN.(FLD), *or* FUNCTION_HANDLE returning an index vector
% from a time series (see, e.g., TS_BOREAL_WARM). Optional ROBUSTACCUMULATOR
% is a FUNCTION_HANDLE used to calculate climatologies (DEFAULT: @nanmean),
% ACCUMULATOR used to compute means and anomalies (DEFAULT: @mean).
%
% Last Saved Time-stamp: <Fri 2010-03-26 17:09:23 Eastern Daylight Time lew.gramer>

  if ( ~exist('fldix','var') || isempty(fldix) )
    fldix = 1:length(stn.(fld).data);
  end;
  if ( isa(fldix,'function_handle') )
    fldix = fldix(stn.(fld));
  end;
  if ( ~exist('robustAccumulator','var') || isempty(robustAccumulator) )
    robustAccumulator = @nanmean;
  end;
  if ( ~exist('accumulator','var') || isempty(accumulator) )
    accumulator = @mean;
  end;
  if ( ~exist('doPlot','var') || isempty(doPlot) )
    doPlot = false;
  end;

  dat = stn.(fld).data(fldix);
  [yr,mo,dy] = datevec(stn.(fld).date(fldix));
  % Floating-point Julian day
  jd = stn.(fld).date(fldix) - datenum(yr,1,1) + 1;

  clims.yrs = unique(yr);
  nyrs = length(clims.yrs);

  % NOTE: We IGNORE data from day 366 of each leap-year
  clims.dtsvec = 1:(1/24):366-(1/24);
  clims.datmtx = repmat(nan, [length(clims.yrs) length(clims.dtsvec)]);

  for yrix = 1:length(clims.yrs)
    ix = find(yr == clims.yrs(yrix));
    [datix,dtsvecix] = intersect_dates(jd(ix),clims.dtsvec,(29.9999/(24*60)));
    % ERROR "vectors must be the same lengths" will occur on the next line,
    % if any date is > 29min 59.64sec from ALL our even one-hour boundaries
    clims.datmtx(yrix,dtsvecix) = dat(ix(datix));
  end;

  fldttl = strrep(fld,'_','\_');

  clims.hrclim = robustAccumulator(clims.datmtx);
  nhrs = numel(clims.hrclim);
  clims.hranom = clims.datmtx - repmat(clims.hrclim, [size(clims.datmtx,1) 1]);
  if ( doPlot )
    figure;
    plot( ((0:nhrs-1)/24)+1, clims.hrclim );
    maxigraph;
    titlename(sprintf('%s hourly %s climatology: %d-%d', ...
                      stn.station_name, fldttl, yr(1), yr(end)));
    xlim([0 366]);
  end;

  ndys = nhrs / 24;
  x = reshape(clims.hrclim, [24 ndys]);
  clims.dyclim = robustAccumulator(x);
  x = reshape(clims.datmtx, [nyrs 24 365]);
  clims.dymean = squeeze(accumulator(x,2));
  clims.dyanom = clims.dymean - repmat(clims.dyclim, [size(clims.dymean,1) 1]);
  if ( doPlot )
    figure;
    plot( 1:ndys, clims.dyclim );
    maxigraph;
    titlename(sprintf('%s daily %s climatology: %d-%d', ...
                      stn.station_name, fldttl, yr(1), yr(end)));
    xlim([0 366]);
  end;

  nwks = floor(ndys / 7);
  x = reshape(clims.dyclim(1:(nwks*7)), [7 nwks]);
  clims.wkclim = robustAccumulator(x);
  x = reshape(clims.dymean(:,1:(nwks*7)), [nyrs 7 nwks]);
  clims.wkmean = squeeze(nanmean(x,2));
  clims.wkanom = clims.wkmean - repmat(clims.wkclim, [size(clims.wkmean,1) 1]);
  if ( doPlot )
    figure;
    plot( 1:nwks, clims.wkclim );
    maxigraph;
    titlename(sprintf('%s weekly %s climatology: %d-%d', ...
                      stn.station_name, fldttl, yr(1), yr(end)));
    xlim([1 52]);
  end;

  nmos = floor(ndys / 30);
  x = reshape(clims.dyclim(1:(nmos*30)), [30 nmos]);
  clims.moclim = robustAccumulator(x);
  clims.moanom = [];
  x = reshape(clims.dymean(:,1:(nmos*30)), [nyrs 30 nmos]);
  clims.momean = squeeze(nanmean(x,2));
  clims.moanom = clims.momean - repmat(clims.moclim, [size(clims.momean,1) 1]);
  if ( doPlot )
    figure;
    plot( 1:nmos, clims.moclim );
    maxigraph;
    titlename(sprintf('%s monthly %s climatology: %d-%d', ...
                      stn.station_name, fldttl, yr(1), yr(end)));
    xlim([1 12]);
  end;

return;
