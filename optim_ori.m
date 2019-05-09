function optim_ori(stn,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX)
%function optim_ori(stn,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX)

  %%%
  %% Call SCRIPT to set:
  %% Variable-name prefixes ("PFX") for various input and output datasets; AND,
  %% All station struct fieldnames used to produce heat budget 
  station_heat_budget_field_names;

  
  switch (lower(stn.station_name))
   case 'fwyf1', minori=-10; maxori= 25;
   case 'mlrf1', minori= 40; maxori= 75;
   case 'mlrf1', minori= 40; maxori= 75;
   case 'lonf1', minori= 60; maxori=100;
   case 'tnrf1', minori= 55; maxori= 95;
   case 'smkf1', minori= 60; maxori=100;
   case 'looe1', minori= 60; maxori=100;
   otherwise,    error('Station coastline orientation: no clue!');
  end;

  if ( isfield(stn,'commentstr') )
    commentstr = stn.commentstr;
  else
    commentstr = '';
  end;

  %%%% DEBUG???
  curInterpMethod = 'linear';
  grdInterpMethod = 'linear';
  commentstr = [commentstr ' ' strrep(KMPFX,'_','\_') ' ' curInterpMethod '/' grdInterpMethod ' '];

  if ( ~isfield(stn,ufld) )
    switch (KMPFX),
      %function stn = get_fkeys_hycom(stn_or_stnm,mindt,maxdt,vars,flds,interpMethod,fkeyspath)
     case 'fkeys_hycom',	stn = get_fkeys_hycom(stn,[],[],[],[],curInterpMethod);
      %function stn = get_gom_hycom(stn_or_stnm,mindt,maxdt,vars,flds,baseurl)
     case 'gom_hycom',	stn = get_gom_hycom(stn);
     otherwise,		error('Unknown km-scale model "%s"',KMPFX);
    end;
    more off;
  end;
  if ( ~isfield(stn,hufld) )
    disp('Recalc hourly');
    stn.(hufld) = interp_ts(stn.(ufld));
    stn.(hvfld) = interp_ts(stn.(vfld));
  end;

  if ( ~isfield(stn,netufld) )
    disp('Recalc net');
    stn.(netufld) = ts_op(stn.(tufld),stn.(hufld),'+');
    stn.(netvfld) = ts_op(stn.(tvfld),stn.(hvfld),'+');
    stn = calc_quasi_eulerian(stn,STOKESPFX,KMPFX,QEPFX);
  end;


  %%%% DEBUG???
  % Run it once with default orientation
                                      % qeufld,qevfld,Tfld,kmtfld,...
                                      % hufld,hvfld,Tfld,kmtfld,...
                                      % ssufld,ssvfld,Tfld,kmtfld,...
  stn = station_cross_shore_advection(stn,bathorifld,...
                                      qeufld,qevfld,Tfld,kmtfld,...
                                      ['raw_' udTfld],udTfld,...
                                      qtfld,qtAdvfld,grdInterpMethod);

  stn = station_calc_kdel2t(stn,model_K_theta,Tfld,...
                            ['raw_' kd2Tfld],kd2Tfld,...
                            qtAdvfld,dTfld,grdInterpMethod);
  stn = station_heat_flux_term_inverse(stn,dTffld,...
                                       dTfld,sfld,[],mhfld);
  stn.(bdTfld) = ts_op(stn.(dTfld),stn.(qbotfld),'+');
  %%%% DEBUG???


  dts = stn.(udTfld).date;
  qdts = stn.(qtAdvfld).date;
  ddts = stn.(dTfld).date;
  bdts = stn.(bdTfld).date;

  [ig,tix] = intersect_dates(dts,stn.(sfld).date);
  tdts = stn.(sfld).date(tix);
  t = stn.(sfld).data(tix);

  [ig,sqtix] = intersect_dates(dts,stn.(sqtfld).date);
  sdts = stn.(sqtfld).date(sqtix);
  sqt = t(1) + cumsum(stn.(sqtfld).data(sqtix));

  oris = minori:5:maxori;
  udT = repmat(nan,[length(dts),length(oris)]);
  qtAdv = repmat(nan,[length(qdts),length(oris)]);
  dT = repmat(nan,[length(ddts),length(oris)]);
  bdT = repmat(nan,[length(bdts),length(oris)]);
  for ix=1:length(oris)
    ori = oris(ix);

    %                                     % qeufld,qevfld,Tfld,kmtfld,...
    %                                     % hufld,hvfld,Tfld,kmtfld,...
    %                                     % ssufld,ssvfld,Tfld,kmtfld,...
    stn = station_cross_shore_advection(stn,ori,...
                                        qeufld,qevfld,Tfld,kmtfld,...
                                        ['raw_' udTfld],udTfld,...
                                        qtfld,qtAdvfld,grdInterpMethod);

    stn = station_calc_kdel2t(stn,model_K_theta,Tfld,...
                              ['raw_' kd2Tfld],kd2Tfld,...
                              qtAdvfld,dTfld,grdInterpMethod);
    stn = station_heat_flux_term_inverse(stn,dTffld,...
                                         dTfld,sfld,[],mhfld);

    udT(:,ix) = stn.(udTfld).data;
    qtAdv(:,ix) = stn.(qtAdvfld).data;
    dT(:,ix) = stn.(dTfld).data;

    stn.(bdTfld) = ts_op(stn.(dTfld),stn.(qbotfld),'+');
    stn = station_heat_flux_term_inverse(stn,bdTffld,bdTfld,sfld,[],mhfld);
    bdT(:,ix) = stn.(bdTfld).data;

  end;


  % res = udT; resnm='\Sigma [ u_k_m ^.\nablaT_k_m ]';
  % res = qtAdv; resnm='\Sigma [ u_k_m ^.\nablaT_k_m + Q_0/\rhoC_ph ]';
  % res = dT; resnm='\Sigma [ u_k_m ^.\nablaT_k_m + K_H_\theta\nabla^2T_k_m + (Q_0)/\rhoC_ph ]';
  res = bdT; resnm='\Sigma [ u_k_m ^.\nablaT_k_m + K_\theta_H\nabla^2T_k_m + (Q_0(\gamma)+Q_b)/\rhoC_ph ]';

  for ix=1:length(oris)
    goodix = find( isfinite(res(:,ix)) );
    res(goodix,ix) = t(1) + cumsum(res(goodix,ix));
  end;

  figure;
  maxigraph;
  hold on;
  plot(tdts,t,'k');
  plot(sdts,sqt,'y');
  % lh = plot(dts,res);
  % lh = plot(qdts,res);
  % lh = plot(ddts,res);
  lh = plot(bdts,res);
  legend(lh, num2str(oris'), 'Location','Best');
  datetick3;
  grid on;
  titlename([commentstr ' ' stn.station_name ': ' resnm]);

return;
