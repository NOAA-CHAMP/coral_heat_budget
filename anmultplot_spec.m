1;

%stnm='smkf1'
%stnm='mlrf1';

for cstnm={'smkf1','lonf1'}
  stnm=cstnm{:};
  stn=[];

  if ( ~exist('stn','var') )
    stn = get_station_from_station_name(stnm);
    stn = load_all_ndbc_data(stn);
  end;
  if ( ~isfield(stn,'Ts') )
    [stn.Ts,stn.Ta,stn.W,stn.P]=intersect_tses(stn.ndbc_sea_t,stn.ndbc_air_t,stn.ndbc_wind1_speed,stn.ndbc_barom);
    stn.nP.date=stn.P.date; stn.nP.data=stn.P.data-1000;

    gapix = [1;find(diff(stn.ndbc_sea_t.date)>(2.1/24))+1;length(stn.ndbc_sea_t.date)];
    [biglen,bigix]=max(stn.ndbc_sea_t.date(gapix(2:end)-1)-stn.ndbc_sea_t.date(gapix(1:end-1)));
    stn.contTs=subset_ts(stn.ndbc_sea_t,gapix(bigix):gapix(bigix+1)-1);
    [ig,stn.contTa,stn.contW,stn.contP]=intersect_tses(stn.contTs,stn.ndbc_air_t,stn.ndbc_wind1_speed,stn.ndbc_barom);

    stn.interpTs = interp_ts(stn.contTs);
    [ig,stn.interpTa,stn.interpW,stn.interpP]=intersect_tses(stn.interpTs,stn.ndbc_air_t,stn.ndbc_wind1_speed,stn.ndbc_barom);
  end;

  doNorm = false;
  %doNorm = true;

  [Pxxes,Wes,fh,lhs]=plot_spec(stn,{'P','W','Ta','Ts'},[],[],[],[],[],doNorm,false);
  set(lhs(end),'LineWidth',1.5);
  set(lhs(end),'Color','r');
  set(lhs(end-1),'Color','k');

  fmg;
  wt_ts(stn.ndbc_sea_t);

end;
