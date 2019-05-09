function [runlen,starti,sgnval] = segrl_ts(ts,tol,dtol)
%function [runlen,starti,sgnval] = segrl_ts(ts,tol,dtol)
%
% Run-length encode the sign of contiguous segments in the time series TS.
% "Sign" here means 0 for (-TOL < TS.data < +TOL), -1 or +1 otherwise, where
% TOL DEFAULTS to 0 (i.e., behavior like function SIGN). "Contiguous segment"
% is defined as one where the gap between successive timestamps is less than
% or equal to DTOL (DEFAULT: 12/24, half a day). In between such gaps, a true
% contiguous time series is ensured with INTERP1. Returned values are RUNLEN,
% vector of integer run-lengths, STARTI, vector of starting indices for each
% run, and SGNVAL, vector of values in {-1,0,+1} corresponding to each run.
%
% A common use for SEGRL_TS would be to segment a heat flux time series into
% periods of cooling, transition, and warming: each element of RUNLEN would
% be the length (in hours, for an hourly time series) of the corresponding
% cooling, transition, or warming period in SGNVAL.
%
% TOL may also be a function handle, accepting a vector of values (TS.data),
% and returning a sequence of indicial values the same length as TS.data. The
% returned values are then a run-length encoding of the return vector of TOL
% (and the set of unique values in vector SGNVAL are those returned by TOL.)
%
% Last Saved Time-stamp: <Sun 2010-07-04 13:36:13 Eastern Daylight Time gramer>


  if ( ~exist('tol','var') || isempty(tol) )
    tol = 0;
  end;
  if ( ~exist('dtol','var') || isempty(dtol) )
    dtol = (12/24);
  end;

  starti = [];
  runlen = [];
  sgnval = [];

  % if ( isa(tol,'function_handle') )
  %   sg = tol(ts.data);
  % else
  %   % Sort of like sg=sign(ts.data) but with tolerance around zero
  %   sg=repmat(0,size(ts.data));
  %   sg(ts.data > +tol) = +1;
  %   sg(ts.data < -tol) = -1;
  % end;

  % Work around big (>dtol) gaps in the time series
  d = diff(ts.date);
  % Goin thru the Big D an I don't mean Dallas
  bigd = [1 find(d > dtol)'];

  for ix = 1:(length(bigd)-1)
    begix = bigd(ix);
    endix = bigd(ix+1)-1;
    dts = ts.date(begix):(1/24):ts.date(endix);
    dat = interp1(ts.date(bigd(ix):bigd(ix+1)),ts.data(bigd(ix):bigd(ix+1)),dts);

    if ( isa(tol,'function_handle') )
      sg = tol(dat);
    else
      % Sort of like sg=sign(ts.data) but with tolerance around zero
      sg = repmat(0, size(dat));
      sg(dat > +tol) = +1;
      sg(dat < -tol) = -1;
    end;

    % (Run-length encoding algorithm courtesy "MATLAB array
    % manipulation tips and tricks", Peter Acklam, Norway)
    kern = find(round(sg(1:end-1)) ~= round(sg(2:end)));
    rl = diff([0 kern(:)' length(sg)]);
    si = bigd(ix) + cumsum(rl);
    ez = sg([kern(:)' length(sg)]);

    runlen(end+1:end+length(rl),1) = rl;
    starti(end+1:end+length(si),1) = si;
    sgnval(end+1:end+length(ez),1) = ez;
  end;

return;
