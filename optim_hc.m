function stn = optim_hc(stn,mos,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX)
%function stn = optim_hc(stn,mos,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX)

 datapath = get_thesis_path('../data');
 figspath = get_thesis_path('../figs');

 %%%
 %% Call SCRIPT to set:
 %% Variable-name prefixes ("PFX") for various input and output datasets; AND,
 %% All station struct fieldnames used to produce heat budget 
 station_heat_budget_field_names;

 opts = station_heat_budget_options(stn);

 for csclg = {'SS','US','SU','UU'}
  sclg = csclg{:};

  %%%
  %% Horizontal Convection

  bet = stn.(slopefld);
  % stn = station_horizontal_convection(stn,sfld,[],mhfld,tspdfld,bq0fld,bet,opts);

  [tix,hix,tspdix,aix,Wix,sq0ix,bq0ix,dTix] = ...
      intersect_all_dates([],stn.(sfld).date,stn.(mhfld).date,stn.(tspdfld).date,stn.(afld).date,stn.(Wfld).date,stn.(sq0fld).date,stn.(bq0fld).date,stn.(bdTffld).date);
  t = stn.(sfld).data(tix);
  s = repmat(36,size(stn.(sfld).data(tix)));
  h = stn.(mhfld).data(hix);
  tspd = stn.(tspdfld).data(tspdix);
  at = stn.(afld).data(aix);
  W = stn.(Wfld).data(Wix);
  sq0 = stn.(sq0fld).data(sq0ix);
  bq0 = stn.(bq0fld).data(bq0ix);
  dT = stn.(bdTffld).data(dTix);

  qstr=bdTffld;
  %%%% ??? DEBUG: Base horizontal convection on surface heating only!
  % q=bq0; stn.commentstr = [stn.commentstr ' HC(Q0+Qb) '];
  %%%% ??? DEBUG: Base horizontal convection on total (km-scale) budget
  q=dT;

  opts.R = get_opt(opts,'R',(1.00-0.08));
  % opts.R = 1.00;

  opts.hc_scaling = sclg;
  % opts.hc_scaling = get_opt(opts,'hc_scaling','SS');
  % opts.hc_scaling = get_opt(opts,'hc_scaling','US');
  % opts.hc_scaling = get_opt(opts,'hc_scaling','SU');
  % opts.hc_scaling = get_opt(opts,'hc_scaling','UU');
  maxhrs = 12;
  opts.maximum_onset_time = get_opt(opts,'maximum_onset_time',maxhrs*3600);
  %%%% ??? DEBUG
  commentstr = [stn.commentstr ' (HC ' opts.hc_scaling ' maxT: ' num2str(maxhrs) 'h) '];
  res = horizontal_convection(t,s,h,q,bet,opts);

  dts = stn.(sfld).date(tix);
  flds = fieldnames(res);
  for fldix = 1:length(flds)
    fld = flds{fldix};
    dat = res.(fld);
    res.(fld) = [];
    res.(fld).date = dts;
    res.(fld).data = dat;

    stnfld = [HCPFX '_' fld];
    stn.(stnfld).date = dts;
    stn.(stnfld).data = dat;
  end;


  % % Simple year-by-year comparison (SUBPLOTs) of truth and heat budget
  % annsubs(stn,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,commentstr,mos);
  % figfname= fullfile(figspath,[lower(stn.station_name) '_' hcdTdt '_annsubs.']);
  % % print('-dtiff',[figfname 'tiff']);
  % % print('-dpng',[figfname 'png']);


  ms1_clim(stn,'daily',RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,commentstr,[],[],mos);
  figfname= fullfile(figspath,[lower(stn.station_name) '_' hcdTdt '_weekly.']);
  % % print('-dtiff',[figfname 'tiff']);
  % % print('-dpng',[figfname 'png']);

 end;

 if ( nargout < 1 )
   stn = []; clear stn;
 end;

return;
