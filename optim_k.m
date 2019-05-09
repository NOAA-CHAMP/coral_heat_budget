function stn = optim_k(stn_or_stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX)
%function stn = optim_k(stn_or_stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX)
%
% Last Saved Time-stamp: <Wed 2011-06-01 18:59:12  lew.gramer>

  grdInterpMethod = 'linear';
  % grdInterpMethod = 'nearest';

  % facs = [0:1:10];
  % facs = [0:.1:1] + 1;
  % facs = [0:.2:1 2:2:10];
  % facs = [0:.2:1 1.5:0.5:4.0];
  facs = [0:10];

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

  stn = station_calc_kdel2t(stn,model_K_theta,Tfld,...
                            ['raw_' kd2Tfld],kd2Tfld,...
                            qtAdvfld,dTfld,grdInterpMethod);

  x.dif = stn.(kd2Tfld);
  x.sfc_btm_adv_dif = ts_op(x.sfc_btm_adv,x.dif,'+');

  % [qix,aix] = intersect_dates(x.sfc_btm.date,x.dif.date);
  [qix,aix] = intersect_dates(x.sfc_btm_adv.date,x.dif.date);
  dts = x.dif.date(aix);
  dat = repmat(nan,[length(dts) length(facs)]);

  for facix = 1:length(facs)
    fac = facs(facix);
    % dat(:,facix) = x.sfc_btm.data(qix) + (fac.*x.dif.data(aix));
    dat(:,facix) = x.sfc_btm_adv.data(qix) + (fac.*x.dif.data(aix));
  end;

  [tix,difix] = intersect_dates(stn.(sfld).date,dts);
  x.t.date = stn.(sfld).date(tix);
  x.t.data = stn.(sfld).data(tix);
  T0 = x.t.data(1);


  [advix,difix] = intersect_dates(x.sfc_btm_adv.date,dts);

  fmg;
  plot(x.t.date,x.t.data,'k-');
  plot(x.sfc_btm_adv.date(advix),T0+cumsum(x.sfc_btm_adv.data(advix)),'g');
  lh = plot(dts,T0+cumsum(dat));
  datetick3;
  legend(lh,num2str(model_K_theta*facs'));
  titlename(strrep(sprintf('%s optimized %s',stn.station_name,kd2Tfld),'_','\_'));


  fmg;
  suptitlename(['Optimized ' upper(stn.station_name) ' ' strrep(sfld,'_','\_') ' vs ' strrep(bdTfld,'_','\_') ' (' grdInterpMethod ')']);
  yrs=2003:2010;
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
  legend(plh,num2str(model_K_theta*facs'),'orientation','hori');

  if ( nargout < 1 )
    stn = []; clear stn;
  end;

return;
