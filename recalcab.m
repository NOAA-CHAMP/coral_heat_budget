function stn = recalcab(stn_or_stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX)
%function stn = recalcab(stn_or_stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX)

  if ( ischar(stn_or_stnm) )
    stnm = stn_or_stnm;
    stn.station_name = stnm;
  else
    stn = stn_or_stnm;
  end;
  if ( ~isfield(stn,'station_name') )
    error('Station struct STN had no .station_name field!');
  end;


  set_more off
  tic,

  datapath = get_thesis_path('../data');
  figspath = get_thesis_path('../figs');

  commentstr = '';

  %%%
  %% Call SCRIPT to set:
  %% Variable-name prefixes ("PFX") for various input and output datasets; AND,
  %% All station struct fieldnames used to produce heat budget 
  station_heat_budget_field_names;

  %%%
  %% Check for presaved results

  matfname = fullfile(datapath,[lower(stn.station_name) '_' TURPFX '_' QEPFX '_heat_budget.mat']);


  if ( ~isfield(stn,dTfld) )
    disp(['Loading from ' matfname]);
    load(matfname,'station');
    flds = fieldnames(station);
    for fldix=1:length(flds)
      fld = flds{fldix};
      stn.(fld) = station.(fld);
    end;
    station = []; clear station;
  end;

%%%% DEBUG??? Fixups: stuff I may have added after the last recalc-from-scratch
    if ( ~isfield(stn,bathorifld) )
      stn = station_optimal_isobath_orientation(stn);
    end;
    if ( ~isfield(stn.(Tfld),'gradient_x') )
      % Calculate gradients and field Laplacians
      stn = calc_field_terms(stn,Tfld);
    end;
    if ( ~isfield(stn,tspdfld) || ~isfield(stn,tdirfld) )
      switch (TIDEPFX),
       case 'tmd_tide',		stn = station_tmd_tide(stn);
       case 'tpxo_tide',	stn = station_tmd_tide(stn);
      end;
    end;
    if ( ~isfield(stn,mhfld) )
      % Calculate the mean tide depth experienced by a watermass moving over
      % an M2 tidal ellipse centered on the coordinates of our station
      stn = station_mean_tide_height(stn,mhfld,bathyfld,hfld);
    end;
%%%% DEBUG???

  parfld = 'bic_surf_par';
  if ( ~isfield(stn,parfld) )
    if ( ~isfield(stn,'sea_t') )
      stn = load_station_data(stn);
    else
      disp('No BIC data!');
    end;
  end;


    %%%
    %% Radiative Fluxes

sai_opts = [];
%%%% ??? DEBUG
% % for Kd=[0.05 0.10 0.20 0.30];
% % for Ab=[0.1 0.2:0.2:1.0];
% %% for Kd=[0.10];
% Kd=[0.30];
% Ab=[];
% % for Ab=[0.1350 0.2350 0.3350]; %0.2350=value calculated for FWYF1,MLRF1,SMKF1
% % for Ab=[0.0690 0.1690 0.2690]; %0.1690=value calculated for LONF1
% commentstr = [' (Kd:' num2str(Kd) ' Ab:' num2str(Ab) ') '];
% sai_opts.kd = Kd;
% sai_opts.bottom_reflectance = Ab;
%%%% ??? DEBUG

      stn = station_absorbed_insolation(stn,asrfld,srfld,mhfld,[],[],gamfld,qbfld,sai_opts);

      % if ( ~isfield(stn,lrfld) )
        % If reanalysis long-wave flux was not specified, user wants a bulk estimate
        %station_bulk_longwave(stn,afld,qfld,pfld,dsrfld,sfld,cfld,dlrf,ulrf,lrf)
        %stn = station_bulk_longwave(stn,afld,qafld,pfld,dsrfld,sfld,cfld,dlrfld,ulrfld,lrfld);
        stn = station_bulk_longwave(stn,afld,qafld,pfld,cfld,sfld,cfld,dlrfld,ulrfld,lrfld);
      % end;

      % Fluxes WITH warm-layer adjustment
      % stn = station_heat_flux(stn,Wfld,afld,rhfld,...
      %                         pfld,sfld,asrfld,lrfld,TURPFX,dsrfld,dlrfld,rfld,...
      %                         Dfld,qeufld,qevfld,wpfld,whfld,pblzfld,true);
      [swix,lwix,lhix,shix,rhix] = ...
          intersect_all_dates([],stn.(asrfld).date,stn.(lrfld).date,...
                              stn.(qlhfld).date,stn.(qshfld).date,stn.(qrhfld).date);
      stn.(q0fld).date = stn.(qlhfld).date(lhix);
      stn.(q0fld).data = stn.(asrfld).data(swix) + stn.(lrfld).data(lwix) ...
          + stn.(qlhfld).data(lhix) + stn.(qshfld).data(shix) + stn.(qrhfld).data(rhix);

      % Algorithm sometimes returns complex numbers!
      stn.(q0fld).data = real(stn.(q0fld).data);
      stn = station_heat_flux_term(stn,q0fld,qtfld,sfld,[],mhfld);

    %%%
    %% Eulerian (km + Stokes) Heat Advection, Km-scale Heat Diffusion

    % stn = station_calc_udotdelt(stn,qeufld,qevfld,Tfld,kmtfld,...
    %                             ['raw_' udTfld],udTfld,...
    %                             qtfld,qtAdvfld);

    % Ignore the along-shore component of model heat advection
    stn = station_cross_shore_advection(stn,bathorifld,...
                                        qeufld,qevfld,Tfld,kmtfld,...
                                        ['raw_' udTfld],udTfld,...
                                        qtfld,qtAdvfld);
    stn = station_calc_kdel2t(stn,K_theta,Tfld,...
                              ['raw_' kd2Tfld],kd2Tfld,...
                              qtAdvfld,dTfld);
    stn = station_heat_flux_term_inverse(stn,dTffld,...
                                         dTfld,sfld,[],mhfld);


    %%%
    %% Benthic Heat Exchanges

    stn = station_benthic_exchange(stn,sfld,netufld,netvfld,qbfld,btfld,qbofld);
    stn = station_heat_flux_term(stn,qbofld,qbotfld,sfld,[],mhfld);

    stn.(bq0fld) = ts_op(stn.(q0fld),stn.(qbofld),'+');
    stn = station_heat_flux_term(stn,bq0fld,bq0tfld,sfld,[],mhfld);

    stn.(bdTfld) = ts_op(stn.(dTfld),stn.(qbotfld),'+');
    stn = station_heat_flux_term_inverse(stn,bdTffld,bdTfld,sfld,[],mhfld);


    %%%
    %% Horizontal Convection

    hc_opts = [];

    bet = stn.(slopefld);
    [tix,hix,tspdix,aix,Wix,q0ix,bq0ix,dTix] = ...
      intersect_all_dates([],stn.(sfld).date,stn.(mhfld).date,stn.(tspdfld).date,stn.(afld).date,stn.(Wfld).date,stn.(q0fld).date,stn.(bq0fld).date,stn.(bdTffld).date);
    t = stn.(sfld).data(tix);
    s = repmat(36,size(stn.(sfld).data(tix)));
    h = stn.(mhfld).data(hix);
    at = stn.(afld).data(aix);
    W = stn.(Wfld).data(Wix);
    q0 = stn.(q0fld).data(q0ix);
    bq0 = stn.(bq0fld).data(bq0ix);
    dT = stn.(bdTffld).data(dTix);

    %%%% ??? DEBUG: Base horizontal convection on surface heating only!
    q=bq0; qstr=bq0fld; commentstr = [commentstr ' HC(Q0+Qb) '];
    %%%% ??? DEBUG: Base horizontal convection on total (km-scale) budget
    % q=dT; qstr=bdTffld;

    hc_opts.R = (1.00-0.08);
    % if ( ~isfield(stn,'ngdc_offshore_slope') )
    %   stn = station_ngdc_offshore_slope(stn);
    % end;

    %scalings = {'SS','US','SU','UU'}
    hc_opts.scaling = 'SS';
    % hc_opts.scaling = 'US';
    % hc_opts.scaling = 'SU';
    % hc_opts.scaling = 'UU';
    hc_opts.maximum_onset_time = 9*3600;
    %%%% ??? DEBUG
    commentstr = [commentstr ' (HC ' hc_opts.scaling ' maxT:9h) '];
    res = horizontal_convection(t,s,h,q,bet,hc_opts);

    dts = stn.(sfld).date(tix);
    flds = fieldnames(res);
    for fldix = 1:length(flds)
      fld = flds{fldix};
      dat = res.(fld);
      res.(fld) = [];
      res.(fld).date = dts;
      res.(fld).data = dat;

      stnfld = ['hc_' fld];
      stn.(stnfld).date = dts;
      stn.(stnfld).data = dat;
    end;

    %%%% ??? DEBUG
    % Truncate start of all time series for display purposes
    % startix = find(dts>=datenum(2004,1,1),1);
    startix = find(dts>=datenum(2006,1,1),1);

    %%%% DEBUG???
    ix = find(~isfinite(q0(startix:end)),1,'last');
    if ( ~isempty(ix) )
      startix = startix + ix;
    end;

    dts = dts(startix:end);
    t = t(startix:end);
    h = h(startix:end);
    at = at(startix:end);
    W = W(startix:end);
    q0 = q0(startix:end);
    bq0 = bq0(startix:end);
    dT = dT(startix:end);
    fac = stn.hc_termFactor.data(startix:end);
    dTdt.date = stn.hc_dTdt.date(startix:end);
    dTdt.data = stn.hc_dTdt.data(startix:end);

    bt = stn.(btfld);
    [ig,startix] = min(abs(bt.date-dts(1)));
    [ig,endix] = min(abs(bt.date-dts(end)));
    bt.date = bt.date(startix:endix);
    bt.data = bt.data(startix:endix);
    %%%% ??? DEBUG

    dsr_lpfld = [asrfld '_24_hour_sum'];
    stn = verify_variable(stn,dsr_lpfld);
    [ig,dsr_lpix] = intersect_dates(dts,stn.(dsr_lpfld).date);

    % lhf_lpfld = [qlhfld '_24_hour_sum'];
    % stn = verify_variable(stn,lhf_lpfld);
    % [ig,lhf_lpix] = intersect_dates(dts,stn.(lhf_lpfld).date);

    udT_lpfld = [udTfld];
    stn = verify_variable(stn,udT_lpfld);
    [ig,udT_lpix] = intersect_dates(dts,stn.(udT_lpfld).date);

    %%%% ??? DEBUG
    T0 = t(1);
    fh=figure; maxigraph; hold on;
    plot(dts,t,dts,at,'k:',bt.date,bt.data,'m--',dts,T0+cumsum(q0.*fac),dts,T0+cumsum(bq0.*fac),dts,T0+cumsum(dT.*fac),dTdt.date,T0+cumsum(dTdt.data));
    % plot(stn.(dsr_lpfld).date(dsr_lpix),(stn.(dsr_lpfld).data(dsr_lpix)./1000),stn.(par_lpfld).date(par_lpix),(stn.(par_lpfld).data(par_lpix).*0.5093./1000),'y:','Color',[.8,.8,0]);
    % plot(stn.(dsr_lpfld).date(dsr_lpix),(stn.(dsr_lpfld).data(dsr_lpix)./1000),stn.(lhf_lpfld).date(lhf_lpix),-(stn.(lhf_lpfld).data(lhf_lpix)./1000),'y:','Color',[.8,.8,0]);
    plot(stn.(dsr_lpfld).date(dsr_lpix),(stn.(dsr_lpfld).data(dsr_lpix)./1000),stn.(udT_lpfld).date(udT_lpix),T0+cumsum(stn.(udT_lpfld).data(udT_lpix)),'y:','Color',[.8,.8,0]);
    datetick3('x',2,'keeplimits');
    legend('T_s','T_a','T_b','Q_0/\rhoC_ph','(Q_0+Q_b)/\rhoC_ph','\partial_tT_k_m',['\partial_tT ' hc_opts.scaling],'\Sigma_1_d Q_S_W/1000','\Sigmau^.\nablaT_k_m', 'Location','Best');
    titlename([ commentstr stn.station_name ' ' strrep(qstr,'_','\_') ]);
    disp(commentstr);
    grid on;
    % axis([datenum(2005,4,1),datenum(2005,5,1),16,28]); datetick3('x',2,'keeplimits');
    axis('tight');
    datetick3('x',2,'keeplimits');
    % ylim([-50 150]);
    % ylim([10 30]);
    ylim([0 50]);
    %%%% ??? DEBUG


    % [stn,hls,hxs]=multiplot_station(stn,{'absorbed_erai_srf','ndbc_erai_30a_latent_heat_flux','ndbc_erai_30a_net_heat_flux','benthic_erai_qbo','benthic_ndbc_erai_30a_net_heat_flux','ndbc_air_t','ndbc_sea_t'},[],[],{'SW','LH','Q_0','BH','Q_0+Q_b','T_a','T_s'},[datenum(2005,4,1),datenum(2005,5,1)],{[-1000,1000],[-1000,1000],[-1000,1000],[-1000,1000],[-1000,1000],[16,28],[16,28]});
    % axes(hxs(1));
    % hold on; plot(stn.erai_srf.date,stn.erai_srf.data,'ro','MarkerSize',1.5);


%%%% ??? DEBUG
% end;
% end;
%%%% ??? DEBUG


  toc,
  set_more;

%%%% ??? DEBUG
close(fh);
chkyr(stn);

  if (nargout < 1)
    stn = [];    clear stn;
  end;

return;
