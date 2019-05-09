function [P,W,Pc] = ch2_spec(ts,N,perDay,CI,plotnm)
%function [P,W,Pc] = ch2_spec(ts,N,perDay,CI,plotnm)
%
% Plot band-averaged, variance-preserving power-spectral density estimate
% for time series TS, averaging over every N frequencies (DEFAULT: 6).
% Normalizes (removes average) from TS first.  If PERDAY is true (DEFAULT),
% display frequencies as 1/day and show super-diurnal through interannual;
% otherwise, show frequencies as 1/hour and focus on near-diurnal/inertial.
% If N is non-finite, averages over a set of frequences of geophysical
% interest - co-tidal, tidal, inertial, weather-band, semi-, and annual; if
% N is a vector, averages over each of those frequencies instead. CI is the
% confidence interval to plot as dashed lines with the PSD (DEFAULT: 0.95).
% PSD plot TITLE and NAME begin with string PLOTNM (DEFAULT: "PMTM").
%
% CALLS: PMTM (Signal Processing Toolbox); NANMEAN (Statistics); PLOT.
%
% Last Saved Time-stamp: <Fri 2013-12-20 17:05:13 Eastern Standard Time gramer>

  if ( ~exist('N','var') || isempty(N) )
    N = 6;
  end;
  if ( ~exist('perDay','var') || isempty(perDay) )
    perDay = true;
  end;
  if ( ~exist('CI','var') || isempty(CI) )
    CI = 0.95;
  end;
  if ( ~exist('plotnm','var') || isempty(plotnm) )
    plotnm = 'PMTM';
  end;

  if ( ~isscalar(N) )
    F = N;
    N = 1;
    [Pxx,Pxxc,W] = pmtm(ts.data-nanmean(ts.data),4,F,1/24,'adapt',CI);
  elseif ( ~isfinite(N) )
    N = 1;
    % "A range of frequencies of geophysical interest..."
    F = (2*pi/24)./([3,6,8,12,12.42,24,24.833,27.4,36,48,24*3,24*6,24*13.66,24*27.32,24*42,24*64,24*182,24*365,3*24*365,6*24*365]/24);
    %F = interp1(1:numel(F),F,1:0.3333333:numel(F))';
    F = interp1(1:numel(F),F,1:0.1000000:numel(F))';
    clear Pxx Pxxc W W_day
    [Pxx,Pxxc,W] = pmtm(ts.data-nanmean(ts.data),4,F,1/24,'adapt',CI);
  else
    %[Pxx,Pxxc,W] = pmtm(ts.data);
    [Pxx,Pxxc,W] = pmtm(ts.data-nanmean(ts.data),[],[],2*pi,'adapt',CI);
  end;

  W_day = (2*pi/24) ./ W;
  W_hour = (2*pi/1) ./ W;

  if ( perDay )
    W_per = W_day;
    % Focus on faster-than-annual periodicities
    %lms = [2e-1,3e3,0,60];
    lms = [2e-2,1e4,0,60];
  else
    % perHour
    W_per = W_hour;
    % Focus on diurnal and inertial periodicities
    lms = [10,40,0,10];
  end;

  % Variance-preserving
  Pxx = Pxx./W_per;
  Pxxc(:,1) = Pxxc(:,1)./W_per;
  Pxxc(:,2) = Pxxc(:,2)./W_per;

  if ( N == 1 )
    W = W_per;
    P = Pxx;
    Pc = Pxxc;
  else
    ix=1:(floor(numel(Pxx)/N)*N);
    midBand = ceil(N/2);
    W = W_per(ix(midBand:N:end));
    P = nanmean(reshape(Pxx(ix),[N,numel(ix)/N]));
    Pc(:,1) = nanmean(reshape(Pxxc(ix,1),[N,numel(ix)/N]));
    Pc(:,2) = nanmean(reshape(Pxxc(ix,2),[N,numel(ix)/N]));
    clear ix;
  end;

  fmg;
  %plot(W,P,'k-',W,Pc(:,1),'k:',W,Pc(:,2),'k:');
  plot(W,P,'k-','LineWidth',2);
  plot(W,Pc(:,1),'k:','LineWidth',3);
  plot(W,Pc(:,2),'k:','LineWidth',3);
  set(gca,'FontSize',16);
  set(gca,'TickLength',[.05,.025]);
  if ( perDay )
    %set(gca,'XScale','log','YScale','log');
    set(gca,'XScale','log','YScale','linear');
    xlabel('Days');
  else
    %set(gca,'XScale','linear','YScale','log');
    set(gca,'XScale','linear','YScale','linear');
    xlabel('Hours');
  end;
  %set(gca,'Position',[0.05,0.10,0.92,0.78]);
  set(gca,'Position',[0.07,0.10,0.90,0.78]);

  axis(lms);
  %grid minor;
  grid on; set(gca,'XMinorGrid','off','YMinorGrid','off');
  titlename([plotnm,', CI=',num2str(roundn(CI*100,-1)),'%, N=',num2str(N)]);

  if ( nargout < 3 )
    Pc=[]; clear Pc;
    if ( nargout < 2 )
      W=[]; clear W;
      if ( nargout < 1 )
        P=[]; clear P;
      end;
    end;
  end;

return;
