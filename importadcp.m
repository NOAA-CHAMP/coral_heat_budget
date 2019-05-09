function [dat,rawdat] = importadcp(fname)
%function [dat,rawdat] = importadcp(fname)
%
% Import an ASCII ADCP datafile FNAME and return the data struct DAT.CURRPROF
% with the timestamped current profiles from the file.
%
% Last Saved Time-stamp: <Mon 2010-07-26 17:45:16 Eastern Daylight Time gramer>

  rawdat = importdata(fname);

  % NaN from IMPORTDATA means that column did not exist on that line!
  % So each non-NaN in rightmost column indicates start of a new profile
  enslns = find(~isnan(rawdat(:,end)));
  if ( enslns(1) ~= 1 )
    error('Bad file format?? No profile header at line 1!');
  end;

  n = length(enslns);
  nlins = size(rawdat,1);
  nbins = nlins - enslns(end);
  % The first NaN in a "profile" row is the first missing column
  endprofln = find(isnan(rawdat(2,:)),1) - 1;

  % Clean up bogus "magic" values
  rawdat(rawdat==999.9) = nan;

  dat.currprof.date = repmat(nan,[n 1]);
  dat.currprof.i_depth = repmat(nan,[n 1]);
  dat.currprof.seatemp = repmat(nan,[n 1]);
  dat.currprof.prof = repmat(nan,[n nbins endprofln]);
  dat.currprof.u = repmat(nan,[n 1]);
  dat.currprof.uprof = repmat(nan,[n nbins]);
  dat.currprof.v = repmat(nan,[n 1]);
  dat.currprof.vprof = repmat(nan,[n nbins]);
  dat.currprof.profgood = repmat(nan,[n nbins]);

  for ix = 1:n
    ensix = enslns(ix);
    if ( ix < n ); endix = enslns(ix+1)-1; else; endix = nlins; end;

    dts = rawdat(ensix,2:7);
    dts(:,1) = 2000 + dts(:,1);
    dat.currprof.date(ix) = datenum(dts);

    dat.currprof.i_depth(ix) = rawdat(ensix,11);
    dat.currprof.seatemp(ix) = rawdat(ensix,12);

    dat.currprof.prof(ix,:,:) = rawdat(ensix+1:endix,1:endprofln);

    dat.currprof.uprof(ix,:) = dat.currprof.prof(ix,:,2) ./ 100;
    dat.currprof.vprof(ix,:) = dat.currprof.prof(ix,:,3) ./ 100;
    dat.currprof.profgood(ix,:) = dat.currprof.prof(ix,:,8);

    % badix = find(dat.currprof.prof(ix,:,8) < 70);
    % dat.currprof.uprof(ix,badix) = nan;
    % dat.currprof.vprof(ix,badix) = nan;

    dat.currprof.u(ix) = nanmean(dat.currprof.uprof(ix,:));
    dat.currprof.v(ix) = nanmean(dat.currprof.vprof(ix,:));
  end;

  if ( nargout < 2 )
    rawdat = []; clear rawdat;
  end;

return;
