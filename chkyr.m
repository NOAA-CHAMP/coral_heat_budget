function x = chkyr(stn,yr,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX)
%function x = chkyr(stn,yr,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX)

  if nargin<2; yr=2005; end;

  station_heat_budget_field_names;

  x.t.date = stn.(sfld).date(get_year(stn.(sfld).date)==yr);
  x.t.data = stn.(sfld).data(get_year(stn.(sfld).date)==yr);
  x = verify_variable(x,'t_24_hour_maximum');
  x = verify_variable(x,'t_60_hour_lowpass');
  x = verify_variable(x,'t_10_day_maximum');

  x.q0.date = stn.(q0fld).date(get_year(stn.(q0fld).date)==yr);
  x.q0.data = stn.(q0fld).data(get_year(stn.(q0fld).date)==yr);

  x.qt.date = stn.(qtfld).date(get_year(stn.(qtfld).date)==yr);
  x.qt.data = x.t.data(1) + cumsum(stn.(qtfld).data(get_year(stn.(qtfld).date)==yr));

  x.bq0.date = stn.(bq0fld).date(get_year(stn.(bq0fld).date)==yr);
  x.bq0.data = stn.(bq0fld).data(get_year(stn.(bq0fld).date)==yr);

  x.bq0t.date = stn.(bq0tfld).date(get_year(stn.(bq0tfld).date)==yr);
  x.bq0t.data = x.t.data(1) + cumsum(stn.(bq0tfld).data(get_year(stn.(bq0tfld).date)==yr));

  x.bdT.date = stn.(bdTfld).date(get_year(stn.(bdTfld).date)==yr);
  x.bdT.data = x.t.data(1) + cumsum(stn.(bdTfld).data(get_year(stn.(bdTfld).date)==yr));

  x.hc_dTdt.date = stn.hc_dTdt.date(get_year(stn.hc_dTdt.date)==yr);
  x.hc_dTdt.data = x.t.data(1) + cumsum(stn.hc_dTdt.data(get_year(stn.hc_dTdt.date)==yr));

  x.sr.date = stn.(srfld).date(get_year(stn.(srfld).date)==yr);
  x.sr.data = stn.(srfld).data(get_year(stn.(srfld).date)==yr);
  x = verify_variable(x,'sr_24_hour_sum');

  x.asr.date = stn.(asrfld).date(get_year(stn.(asrfld).date)==yr);
  x.asr.data = stn.(asrfld).data(get_year(stn.(asrfld).date)==yr);
  x = verify_variable(x,'asr_24_hour_sum');

  x.lr.date = stn.(lrfld).date(get_year(stn.(lrfld).date)==yr);
  x.lr.data = stn.(lrfld).data(get_year(stn.(lrfld).date)==yr);
  x = verify_variable(x,'lr_24_hour_sum');

  x.qlh.date = stn.(qlhfld).date(get_year(stn.(qlhfld).date)==yr);
  x.qlh.data = stn.(qlhfld).data(get_year(stn.(qlhfld).date)==yr);
  x = verify_variable(x,'qlh_24_hour_sum');
  % x.qlh_24_hour_sum.data = x.qlh_24_hour_sum.data*9e-5;

  x.qsh.date = stn.(qshfld).date(get_year(stn.(qshfld).date)==yr);
  x.qsh.data = stn.(qshfld).data(get_year(stn.(qshfld).date)==yr);
  x = verify_variable(x,'qsh_24_hour_sum');

  x.qrh.date = stn.(qrhfld).date(get_year(stn.(qrhfld).date)==yr);
  x.qrh.data = stn.(qrhfld).data(get_year(stn.(qrhfld).date)==yr);

  x.radif = ts_op(x.qlh,x.qsh,'+');
  x = verify_variable(x,'radif_24_hour_sum');

  x.coolf = ts_op(x.lr,x.radif,'+');
  x = verify_variable(x,'coolf_24_hour_sum');

  x.non_qlh_qb = ts_op(x.q0,x.qlh,'-');
  x.non_qlh = ts_op(x.bq0,x.qlh,'-');

  figure;
  maxigraph;
  grid on;
  ax(1)=subplot_tight(3,1,[1]); hold on;
  plot_ts(x.t,x.bq0t,x.hc_dTdt); ylim([10,40]);
  % plot_ts(x.t_10_day_maximum,x.bq0t); ylim([10,40]);
  % plot_ts(x.t_60_hour_lowpass,x.qt); ylim([10,40]);
  ax(2)=subplot_tight(3,1,[2,3]); hold on;
  % plot_ts(x.asr_24_hour_sum,x.lr_24_hour_sum,x.qlh_24_hour_sum,x.qsh_24_hour_sum,x.coolf_24_hour_sum);
  plot_ts(x.non_qlh,x.qlh); legend(ax(2),'Q_0+Q_b-Q_L_H','Q_L_H');
  % plot_ts(x.non_qlh_qb,x.qlh); legend(ax(2),'Q_0-Q_L_H','Q_L_H');

  if ( nargout < 1 )
    x = []; clear x;
  end;

return;
