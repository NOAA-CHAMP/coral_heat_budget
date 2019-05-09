function stn = ms1_clim(stn,per,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,substitute_field_names,commentstr,tempfld,budgetfld,mos,begyr,endyr)
%function stn = ms1_clim(stn,per,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,substitute_field_names,commentstr,tempfld,budgetfld,mos,begyr,endyr)
%
% Plot comparison of ocean heat budget climatology (expressed as an implied
% sea teperature time series) with climatology of actual sea temperature.
%
% Last Saved Time-stamp: <Tue 2012-07-31 16:56:01  lew.gramer>


  if ( ~exist('per','var') || isempty(per) )
    per = 'weekly';
  end;
  % Allowing leap-days might cause problems when comparing cum stats on
  % *different time series*, e.g., a time series which includes data for
  % one or more 29th's of Feb, and one which does not. If caller does not
  % like this behavior, they can simply specify custom CUMFUN and MINN.
  switch ( per ),
   case 'hourly',   n = 365*24; minN=1;     perstr='_h_r'; N = 1;
   case 'daily',    n = 365;    minN=23;    perstr='_d_y'; N = 24;
   case 'pentad',   n = 73;     minN=23*5;  perstr='_5_d'; N = 24*5;
   case 'weekly',   n = 52;     minN=23*7;  perstr='_w_k'; N = 24*7;
   case 'monthly',  n = 12;     minN=23*28; perstr='_m_o'; N = 24*[31,28,31,30,31,30,31,31,30,31,30,31]';
   case 'seasonal', n = 4;      minN=23*90; perstr='_s_e_a_s'; N = 24*[31+28+31,30+31+30,31+31+30,31+30+31]';
   otherwise,       error('Do not know how to do a "%s"-period climatology!',char(per));
  end;

  if ( ~exist('commentstr','var') || isempty(commentstr) )
    if ( isfield(stn,'commentstr') )
      commentstr = stn.commentstr;
    else
      commentstr = '';
    end;
  end;

  x.commentstr = commentstr;
  station_heat_budget_field_names;
  %% HACK ALERT: string COMMENTSTR is *not* a variable name...
  commentstr = x.commentstr;
  clear x

  if ( ~exist('tempfld','var') || isempty(tempfld) )
    tempfld = sfld;
  end;
  if ( ~exist('budgetfld','var') || isempty(budgetfld) )
    budgetfld = hcdTdt;
  end;
  if ( ~exist('mos','var') || isempty(mos) )
    mos = 1:12;
  end;
  if ( ~exist('begyr','var') || isempty(begyr) )
    begyr = 1900;
  end;
  if ( ~exist('endyr','var') || isempty(endyr) )
    endyr = 2100;
  end;

  % Only applies for PER == 'daily'
  begjd = datenum(1,mos(1),1) - datenum(1,1,1) + 1;
  endjd = datenum(1,mos(end)+1,1) - datenum(1,1,1);


  % Would prefer MEDIAN, but published climatologies generally use MEAN (and
  % for, e.g., daily insolation, these give very different results!)
  sumfun = @nanmean;

  %SEE function [cum,tid] = grp_ts(dat,dts,per_or_cumfun,sumfun,minN)

  begdt = max([datenum(begyr,min(mos(:)),1), ...
               stn.(tempfld).date(find(isfinite(stn.(tempfld).data),1)), ...
               stn.(budgetfld).date(find(isfinite(stn.(budgetfld).data),1)) ]);
  enddt = min([datenum(endyr,max(mos(:))+1,1), ...
               stn.(tempfld).date(find(isfinite(stn.(tempfld).data),1,'last')), ...
               stn.(budgetfld).date(find(isfinite(stn.(budgetfld).data),1,'last')) ]);
  tix = find(begdt<=stn.(tempfld).date & stn.(tempfld).date<=enddt);
  bix = find(begdt<=stn.(budgetfld).date & stn.(budgetfld).date<=enddt);

  x.t = grp_ts(stn.(tempfld).data(tix),stn.(tempfld).date(tix),per,sumfun,minN);
  x.hc_dTdt = N.*cumsum(grp_ts(stn.(budgetfld).data(bix),stn.(budgetfld).date(bix),per,sumfun,minN));

  if ( n == 365 )
    begn = begjd;
    endn = endjd;
  else
    begn = 1;
    endn = n;
  end;

  fmg;
  plot([begn:endn],x.t(begn:endn),'k-',...
       [begn:endn],x.t(begn)+x.hc_dTdt(begn:endn)-x.hc_dTdt(begn),'r-');
  plot([begn:4:endn],x.t(begn:4:endn),'k.',...
       [begn:4:endn],x.t(begn)+x.hc_dTdt(begn:4:endn)-x.hc_dTdt(begn),'ro');
  plh=plot([begn],x.t(begn),'k.-',...
           [begn],x.t(begn)+x.hc_dTdt(begn)-x.hc_dTdt(begn),'ro-');
  legend(plh,'T_s',['T_s(0)+\Sigma' perstr '\partial_tT_H_C'],...
         'Location','SouthEast', 'Orientation','horizontal');
  axis([begn,endn,15,33]);
  yrstr = [ num2str(get_year(begdt)) '-' num2str(get_year(enddt)) ];
  mostr = [ num2str(mos(1)) '-' num2str(mos(end)) ];
  titlename([upper(stn.station_name) ': ' ...
             strrep(budgetfld,'_','\_') ' ' upper(per) ' climatology ' ...
             yrstr ' ' mostr ' ' commentstr]);

  if ( nargout < 1 )
    stn = []; clear stn;
  end;

return;
