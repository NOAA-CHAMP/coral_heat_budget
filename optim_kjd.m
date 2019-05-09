function stn = optim_kjd(stn_or_stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX)
%function stn = optim_kjd(stn_or_stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX)
%
% Last Saved Time-stamp: <Mon 2011-04-18 15:46:52  Lew.Gramer>

  grdInterpMethod = 'linear';
  % peakjds = [0:60:300];
  peakjds = [30 60 90  270 300 330];

  station_heat_budget_field_names;

  if ( ~isfield(stn_or_stnm,bdTffld) )
    stn = station_heat_budget(stn_or_stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX);
  else
    stn = stn_or_stnm;
  end;
  clear stn_or_stnm;

  disp(udTfld);

  x.sfc_btm = stn.(bq0tfld);
  x.sfc_btm_adv = ts_op(x.sfc_btm,stn.(udTfld),'+');

  stn = station_calc_kdel2t(stn,K_theta,Tfld,...
                            ['raw_' kd2Tfld],kd2Tfld,...
                            qtAdvfld,dTfld,grdInterpMethod);

  x.dif = stn.(kd2Tfld);
  x.sfc_btm_adv_dif = ts_op(x.sfc_btm_adv,x.dif,'+');

  % [qix,aix] = intersect_dates(x.sfc_btm_adv.date,x.dif.date);
  [qix,aix] = intersect_dates(x.sfc_btm.date,x.dif.date);
  dts = x.dif.date(aix);
  dat = repmat(nan,[length(dts) length(peakjds)]);

  [tix,difix] = intersect_dates(stn.(sfld).date,dts);
  x.t.date = stn.(sfld).date(tix);
  x.t.data = stn.(sfld).data(tix);
  T0 = x.t.data(1);

  mink = 0.5*K_theta;
  maxk = 7.5*K_theta;
  for jdix = 1:length(peakjds)
    peakjd = peakjds(jdix);
    stn = station_calc_kdel2t(stn,[mink,maxk,peakjd],Tfld,...
                              ['raw_' kd2Tfld],kd2Tfld,...
                              qtAdvfld,dTfld,grdInterpMethod);
    % dat(:,jdix) = x.sfc_btm_adv.data(qix) + stn.(kd2Tfld).data(aix);
    dat(:,jdix) = x.sfc_btm.data(qix) + stn.(kd2Tfld).data(aix);
    %DEBUG:    fmg; plot(dts,T0+cumsum(stn.(kd2Tfld).data(aix))); titlename(num2str(peakjd));
  end;


  fmg;
  plot(x.t.date,x.t.data,'k-');
  plot(x.sfc_btm.date,T0+cumsum(x.sfc_btm.data),'g');
  lh = plot(dts,T0+cumsum(dat));
  datetick3;
  legend(lh,num2str(peakjds'));
  titlename(strrep(sprintf('%s optimized %s',stn.station_name,kd2Tfld),'_','\_'));


  fmg;
  suptitlename([upper(stn.station_name) ' ' strrep(sfld,'_','\_') ' vs ' strrep(bdTfld,'_','\_') ' (' grdInterpMethod ')']);
  yrs=2004:2010;
  for yrix=1:length(yrs)
    yr = yrs(yrix);
    t.date = x.t.date(get_year(x.t.date)==yr);
    t.data = x.t.data(get_year(x.t.date)==yr);

    q.date = dts(get_year(dts)==yr);
    q.data = dat(get_year(dts)==yr,:);
    if ( ~isempty(q.data(isfinite(q.data))) )
      qs.date = q.date;
      qs.data = t.data(1) + cumsum(q.data) - q.data(1);

      subplot_tight(4,2,yrix);
      hold on;
      % Plot real sea temperature for comparison
      plot(t.date,t.data,'k');
      % Plot all experiment resultants at once
      plh = plot(qs.date,qs.data);
      xlim([datenum(yr,1,1),datenum(yr+1,1,1)]);
      datetick('x',2,'keeplimits','keepticks');
      set_datetick_cursor;
      ylim([0,45]);
    end;
  end;
  legend(plh,num2str(peakjds'),'orientation','hori');

  if ( nargout < 1 )
    stn = []; clear stn;
  end;

return;
