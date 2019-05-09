function stn = optim_Kt(stn_or_stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX)
%function stn = optim_Kt(stn_or_stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX)
%
% Last Saved Time-stamp: <Fri 2011-04-15 10:50:50  Lew.Gramer>

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

  grdInterpMethod = 'linear';

  stn = station_calc_kdel2t(stn,K_theta,Tfld,...
                            ['raw_' kd2Tfld],kd2Tfld,...
                            qtAdvfld,dTfld,grdInterpMethod);


  x.dif = stn.(kd2Tfld);
  x.sfc_btm_adv_dif = ts_op(x.sfc_btm_adv,x.dif,'+');

  K_thetas = 0:2.5:25;

  [qix,aix] = intersect_dates(x.sfc_btm_adv.date,x.dif.date);
  dts = x.sfc_btm_adv.date(qix);
  dat = repmat(nan,[length(dts) length(K_thetas)]);

  for kix = 1:length(K_thetas)
    kth = K_thetas(kix);
    stn = station_calc_kdel2t(stn,kth,Tfld,...
                              ['raw_' kd2Tfld],kd2Tfld,...
                              qtAdvfld,dTfld,grdInterpMethod);
    dat(:,kix) = x.sfc_btm_dif.data(qix) + (fac.*x.adv.data(aix));
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
