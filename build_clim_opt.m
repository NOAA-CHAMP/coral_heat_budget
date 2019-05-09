function val = build_clim_opt(op,opnm,dts,doDebug)
%function val = build_clim_opt(op,opnm,dts,doDebug)
%
% Build a term from the option value OP, of suitable length for combination
% with time series matching the vector of timestamps (DATENUMs) DTS. OP may
% either be a time series (with .date and .data fields), a scalar constant,
% a numerical 3-vector to specify an ANNUAL cosine-wave climatology (i.e.,
% OP=[MIN-ANNUAL-VAL,MAX-ANNUAL-VAL,PEAK-YEAR-DAY]), a 4-vector for a sub-
% or super-annual cosine climatology (fourth element must be PERIOD of cosine
% wave in days), or an (Nx4)-vector - to form a climatology from the sum of N
% cosine waves. OP may also be cell array specifying functional relationship
% between the term and some other data source: in this case, the first cell
% must be a time series (struct with .date,.data), the second a function that
% accepts and returns a numeric vector. OPNM is an optional name string for
% option OP, used only for diagnostic and error messages, and in plots.
%
% If optional DODEBUG is not specified or TRUE, diagnostic info is DISP'd.
%
% Last Saved Time-stamp: <Sat 2012-03-31 01:50:53  Lew.Gramer>

  if ( ~exist('doDebug','var') || isempty(doDebug) )
    doDebug = true;
  end;

  if ( isempty(opnm) )
    opnm = inputname(1);
    if ( isempty(opnm) )
      opnm = 'OP';
    end;
  end;
  opnm(1) = upper(opnm(1));

  if ( is_valid_ts(op) && isnumeric(op.data) )

    if ( doDebug )
      disp( [ 'Interpolating time series ' opnm ] );
    end;
    val = interp1(op.date,op.data,dts);
    val(~isfinite(val)) = [];
    if ( numel(val) ~= numel(dts) )
      error('Cannot extrapolate! Date range mismatch between %s time series and DTS',opnm);
    end;

  elseif ( isnumeric(op) )

    switch ( numel(op) )
     case 1,
      % Constant scalar value
      val = real(op);
      if ( doDebug )
        disp( [ 'Using constant ' opnm ' = ' num2str(val) ] );
      end;

     case 3,
      % Seasonally (annually) varying value
      minval   = op(1);
      maxval   = op(2);
      peakjd   = op(3);
      meanval  = mean([minval maxval]);
      valrng   = maxval - meanval;
      peakrads = (2*pi*peakjd/366);
      if ( doDebug )
        disp( [ 'Using annual cosine-modulated ' opnm ' = ' num2str([minval,maxval,peakjd],'%g ') ] );
      end;
      val = meanval + valrng.*cos((2.*pi.*get_yearday(dts)/366)-peakrads);
      %DEBUG:      fmg; plot(dts,val); datetick3;
      %DEBUG:      fmg; [cum,tid]=grp_ts(val,dts,[],[],1); plot(tid,cum); titlename([strrep(opnm,'_','\_'),': ',num2str(op,'%g ')]);

     case 4,
      % Sub- or super-annual climatology (i.e., cosine curve not exactly 1 year long)
      minval   = op(1);
      maxval   = op(2);
      peakjd   = op(3);
      % Period (days) of cyclic variation: May be >366 for INTERANNUAL
      period   = op(4);
      meanval  = mean([minval maxval]);
      valrng   = maxval - meanval;
      peakrads = (2*pi*peakjd/period);
      if ( doDebug )
        disp( [ 'Using non-annual cosine-modulated ' opnm ' = ' num2str([minval,maxval,peakjd,period],'%g ') ] );
      end;
      val = meanval + valrng.*cos((2.*pi.*get_yearday(dts)./period)-peakrads);
      %DEBUG:      fmg; plot(dts,val); datetick3;
      %DEBUG:      fmg; [cum,tid]=grp_ts(val,dts,[],[],1); plot(tid,cum); titlename([opnm,': ',num2str(op,'%g ')]);

     otherwise,
      if ( mod(numel(op),4) ~= 0 )
        error('Numeric %s may be scalar, annual 3-vector, non-annual 4-vector, or (Nx4)-vector!',opnm);
      else
        if ( doDebug )
          disp( [ 'Using multi-cosine-modulated ' opnm ' = ' num2str(op,'%g ') ] );
        end;
        val = repmat(0,size(dts));
        for opix = 1:4:numel(op)
          minval   = op(1);
          maxval   = op(2);
          peakjd   = op(3);
          period   = op(4);
          meanval  = mean([minval maxval]);
          valrng   = maxval - meanval;
          peakrads = (2*pi*peakjd/period);
          val = val + meanval + valrng.*cos((2.*pi.*get_yearday(dts)./period)-peakrads);
        end;
        %DEBUG:        fmg; plot(dts,val); datetick3;
        %DEBUG:        fmg; [cum,tid]=grp_ts(val,dts,[],[],1); plot(tid,cum); titlename([opnm,': ',num2str(op,'%g ')]);
      end;
    end;

  elseif ( iscell(op) && numel(op) == 2 )
    ts = op{1};
    fn = op{2};
    if ( doDebug )
      disp( [ 'Using functional ',char(fn),'(TS)' ] );
    end;
    % [ig,ix] = intersect_dates(dts,ts.date);
    % val = fn(ts.data(ix));
    % NOTE: Have to interpolate regardless of big gaps, so vector lengths match
    dat = interp1(ts.date,ts.data,dts);
    val = fn(dat);
    %DEBUG:    fmg; plot(dts,val); datetick3;
    %DEBUG:    fmg; [cum,tid]=grp_ts(val,dts,[],[],1); plot(tid,cum); titlename([opnm,': ',char(fn)]);

  else
    error('%s must be a valid time series, scalar, 3- or 4-vector, or cell 2-array!',opnm);
  end;

return;
