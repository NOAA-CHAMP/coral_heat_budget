function stn = optim_km(stn_or_stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX)
%function stn = optim_km(stn_or_stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX)
%
% Last Saved Time-stamp: <Tue 2011-04-19 10:39:50  Lew.Gramer>

  station_heat_budget_field_names;

  if ( ~isfield(stn_or_stnm,bdTffld) )
    stn = station_heat_budget(stn_or_stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX);
  else
    stn = stn_or_stnm;
  end;
  clear stn_or_stnm;

  disp(kd2Tfld);

  x.sfc_btm = stn.(bq0tfld);
  x.sfc_btm_dif = ts_op(x.sfc_btm,stn.(kd2Tfld),'+');

  grdInterpMethod = 'linear';
  % stn = station_calc_udotdelt(stn,qeufld,qevfld,Tfld,kmtfld,...
  %                             ['raw_' udTfld],udTfld,...
  %                             [],[],grdInterpMethod);
  stn = station_cross_shore_advection(stn,bathorifld,...
                                      qeufld,qevfld,Tfld,kmtfld,...
                                      ['raw_' udTfld],udTfld,...
                                      [],[],grdInterpMethod);

  x.adv = stn.(udTfld);
  x.sfc_btm_dif_adv = ts_op(x.sfc_btm_dif,x.adv,'+');

  facs = 0:.1:1;

  [qix,aix] = intersect_dates(x.sfc_btm_dif.date,x.adv.date);
  dts = x.sfc_btm_dif.date(qix);
  dat = repmat(nan,[length(dts) length(facs)]);

  for facix = 1:length(facs)
    fac = facs(facix);
    dat(:,facix) = x.sfc_btm_dif.data(qix) + (fac.*x.adv.data(aix));
  end;

  [tix,advix] = intersect_dates(stn.(sfld).date,dts);
  x.t.date = stn.(sfld).date(tix);
  x.t.data = stn.(sfld).data(tix);
  T0 = x.t.data(1);

  fmg;
  plot(x.t.date,x.t.data,'k-');
  plot(x.sfc_btm.date,T0+cumsum(x.sfc_btm.data),'g');
  lh = plot(dts,T0+cumsum(dat));
  datetick3;
  legend(lh,num2str(facs'));
  titlename(strrep(sprintf('%s optimized %s',stn.station_name,udTfld),'_','\_'));

  if ( nargout < 1 )
    stn = []; clear stn;
  end;

return;
