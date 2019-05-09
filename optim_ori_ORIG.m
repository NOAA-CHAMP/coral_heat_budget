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
   case 'lonf1', minori= 60; maxori=100;
   case 'tnrf1', minori= 55; maxori= 95;
   case 'smkf1', minori= 60; maxori=100;
   case 'looe1', minori= 60; maxori=100;
   otherwise,    error('Station coastline orientation: no clue!');
  end;

  dts = stn.(udTfld).date;
  qdts = stn.(qtAdvfld).date;
  bdts = stn.(bdTfld).date;

  oris = minori:5:maxori;
  udT = repmat(nan,[length(dts),length(oris)]);
  qtAdv = repmat(nan,[length(qdts),length(oris)]);
  bdT = repmat(nan,[length(bdts),length(oris)]);
  for ix=1:length(oris)
    ori = oris(ix);

    stn.qeu.date = stn.(qeufld).date;
    [stn.qeu.data,ig] = reorient_vectors(ori,stn.(qeufld).data,stn.(qevfld).data);

    stn.qezero.date = stn.(qeufld).date;
    stn.qezero.data = repmat(0,size(stn.qezero.date));

    stn.Tfld.date = stn.(Tfld).date;
    midx = round(length(stn.(Tfld).lon) / 2);
    midy = round(length(stn.(Tfld).lat) / 2);
    stn.Tfld.lon = stn.(Tfld).lon(midx);
    stn.Tfld.lat = stn.(Tfld).lat(midy);
    stn.Tfld.field = stn.(Tfld).field(:,midx,midy);
    gx = stn.(Tfld).gradient_x(:,midx,midy);
    gy = stn.(Tfld).gradient_y(:,midx,midy);
    [stn.Tfld.gradient_x,ig] = reorient_vectors(ori,gx,gy);
    stn.Tfld.gradient_y = repmat(0,size(stn.Tfld.gradient_x));

    stn.qt = stn.(qtfld);
    stn = station_calc_udotdelt(stn,'qeu','qezero','Tfld',...
                                ['raw_' 'udT'],'udT',...
                                'qt','qtAdv');
    udT(:,ix) = stn.udT.data;
    qtAdv(:,ix) = stn.qtAdv.data;

    stn.bdT = ts_op(stn.qtAdv,stn.(qbotfld),'+');
    stn = station_heat_flux_term_inverse(stn,'bdTf','bdT',sfld,[],mhfld);
    bdT(:,ix) = stn.bdT.data;
  end;


  res = udT; resnm='\Sigma u^.\nablaT';
  res = qtAdv; resnm='\Sigma u^.\nablaT + Q_0/\rhoC_ph';
  res = bdT; resnm='\Sigma u^.\nablaT + (Q_0+Q_b)/\rhoC_ph';

  for ix=1:length(oris)
    goodix = find( isfinite(res(:,ix)) );
    res(goodix,ix) = cumsum(res(goodix,ix));
  end;

  figure;
  maxigraph;
  % plot(dts,res);
  % plot(qdts,res);
  plot(bdts,res);
  legend(num2str(oris'), 'Location','Best');
  datetick3;
  grid on;
  titlename([stn.station_name ': Accumulated heat advection ' resnm]);

return;
