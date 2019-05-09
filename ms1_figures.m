function ms1_figures(stn)
%function ms1_figures(stn)

  stnm = upper(stn.station_name);

  flds = {'ndbc_ncep_30a_heat_flux_term','gom_hycom_dt','netqf','gom_hycom_dt_netqf','gom_hycom_netqf'};
  flds = {'netqf','gom_hycom_netqf'};
  legs = {'HC(Q_0) / week (thermal siphon)', 'HC(Q_0 + GoM 4km HYCOM) / week', };

  for fldix = 1:length(flds)
    fld = flds{fldix};
    %DEBUG:
    disp(fld);

    [yr,mo,dy,hr,mi,sc] = datevec(stn.(fld).date);
    wk = get_week(stn.(fld).date);
    jd = get_jday(stn.(fld).date);
    % [mn,se,n,st,ci] = grpstats(stn.(fld).data,hr,{'mean','sem','numel','std','meanci'});
    % [mn,n,ci] = grpstats(real(stn.(fld).data),mo,{'mean','numel','meanci'});
    [mn,n,ci] = grpstats(real(stn.(fld).data),wk,{'mean','numel','meanci'});
    % [mn,n,ci] = grpstats(real(stn.(fld).data),jd,{'mean','numel','meanci'});
    % mnci = bootweek(stn,fld,30);
    % mnci = bootmon(stn,fld,30);
    % mn = mnci(2,:);
    % if ( size(mnci,1) > 2 )
    %   ci = mnci([1 3],:);
    % end;

    % errorbar(1:length(mn),mn,ci(:,2)-mn);

    % figure;  bar(1:length(mn),n./max(n(:)));

    % mn(n < (0.5*max(n(:)))) = nan;

    % Heating terms are [K/hour] - scale appropriately for comparison below
    switch (length(mn))
     case 12,
      mns(fldix,1:length(mn)) = mn(:)' .* (eomday(2009,[1:12]).*24);
     case 52,
      mns(fldix,1:length(mn)) = mn(:)' .* (7*24);
     case {365,366},
      mns(fldix,1:length(mn)) = mn(:)' .* 24;
     case {365*24,366*24},
      mns(fldix,1:length(mn)) = mn(:)';
     case 24,
      mns(fldix,1:length(mn)) = mn(:)';
     otherwise,
      error('Unhandled sampling period %d',length(mn));
    end;
  end;

  %DEBUG:  figure;  bar(1:length(mn),n./max(n(:)));

  fld = 'ndbc_sea_t';
  % % mnci = bootweek(stn,fld,30);
  % mnci = bootmon(stn,fld,30);
  % mn = mnci(2,:);

  [yr,mo,dy,hr,mi,sc] = datevec(stn.(fld).date);
  wk = get_week(stn.(fld).date);
  jd = get_jday(stn.(fld).date);
  % [mn,n,ci] = grpstats(real(stn.(fld).data),mo,{'mean','numel','meanci'});
  [mn,n,ci] = grpstats(real(stn.(fld).data),wk,{'mean','numel','meanci'});
  % [mn,n,ci] = grpstats(real(stn.(fld).data),jd,{'mean','numel','meanci'});
  flds{end+1} = ['\delta(' fld ')'];
  legs{end+1} = ['\Delta(T) / week'];
  %DEBUG:
  disp(flds{end});
  mns(end+1,1:size(mns,2)) = [ (mn(1)-mn(end)) diff(mn(:))' ];

  figure;
  maxigraph;
  plot(1:size(mns,2),mns);
  % legend(strrep(flds,'_','\_'), 'Location','Best');
  legend(legs, 'Location','Best');
  xlim([1 size(mns,2)]);
  ylabel('^oC / week');
  % print('-dtiff','../figs/fwyf1-weekly-climatology-comparison-netqf-fluxes.tiff');


  flds = {}; mns = [];
  fldbasenms = {'ncep_srf','ncep_lrf', ...
                'monthly_nocs_srf','monthly_nocs_lrf', ...
                'landy_sr','landy_lr'};
  legs = {'NCEP NARR Q_S_W',       'NCEP NARR Q_L_W', ...
          'NOC Southampton Q_S_W', 'NOC Southampton Q_L_W', ...
          'Large&Yeager (2009) Q_S_W', 'Large & Yeager (2009) Q_L_W', ...
         };
  for ix = 1:length(fldbasenms)
    flds{ix} = [fldbasenms{ix}];
    mns(ix,1:12) = grpstats(real(stn.(flds{ix}).data),get_month(stn.(flds{ix}).date));
    mns(ix,1:12) = grpstats(real(stn.(flds{ix}).data),get_month(stn.(flds{ix}).date));
  end;

  figure; maxigraph; hold on;
  plot(1:12,mns); xlim([1 12]);
  % legend(strrep(flds,'_','\_'), 'Location','Best');
  legend(legs, 'Location','Best');
  titlename([stnm ' Radiative Fluxes (NCEP NARR)']);
  ylabel('Monthly Mean Flux [ W / m^2 ]');
  % print('-dtiff','../figs/fwyf1-monthly-climatology-comparison-radiative-fluxes.tiff');

  flds = {}; mns = [];
  fldbasenms = { 'ndbc_ncep_30a_sensible_heat_flux','ndbc_ncep_30a_latent_heat_flux', ...
                 'monthly_nocs_sensible_heat_flux','monthly_nocs_latent_heat_flux', ...
                 'landy_sh','landy_lh' };
  legs = {      'TOGA-COARE/NCEP NARR Q_S_H',    'TOGA-COARE/NCEP NARR Q_L_H', ...
                'NOC Southampton Q_S_H',         'NOC Southampton Q_L_H', ...
                'Large&Yeager (2009) Q_S_H',     'Large & Yeager (2009) Q_L_W', ...
         };
  for ix = 1:length(fldbasenms)
    flds{ix} = [fldbasenms{ix}];
    mns(ix,1:12) = grpstats(real(stn.(flds{ix}).data),get_month(stn.(flds{ix}).date));
  end;
  figure; maxigraph; hold on;
  plot(1:12,mns); xlim([1 12]);
  % legend(strrep(flds,'_','\_'), 'Location','Best');
  legend(legs, 'Location','Best');
  titlename([stnm ' Turbulent Fluxes (TOGA-COARE 3.0a)']);
  ylabel('Monthly Mean Flux [ W / m^2 ]');
  % print('-dtiff','../figs/fwyf1-monthly-climatology-comparison-turbulent-fluxes.tiff');

  flds = {}; mns = [];
  fldbasenms = {'gom_hycom_advected_heat','qf','gom_hycom_qf'};
  legs =       { 'GoM 4km HYCOM u ^. \nabla T', 'HC(Q_0) u ^. \nabla T', ...
                 'HC(Q_0 + GoM 4km HYCOM) u ^. \nabla T'};
  for ix = 1:length(fldbasenms)
    flds{ix} = [fldbasenms{ix}];
    mns(ix,1:12) = grpstats(real(stn.(flds{ix}).data),get_month(stn.(flds{ix}).date));
  end;
  figure; maxigraph; hold on;
  plot(1:12,mns); xlim([1 12]);
  % legend(strrep(flds,'_','\_'), 'Location','Best');
  legend(legs, 'Location','Best');
  titlename([stnm ' Heat Advection (GoM HYCOM, thermal siphon)']);
  ylabel('Advected Heat [ ^oC / hr ]');
  % print('-dtiff','../figs/fwyf1-monthly-climatology-comparison-advected-heat.tiff');


  flds = {}; mns = [];
  fldbasenms = {'ndbc_ncep_30a_net_heat_flux', 'ndbc_hfbulk_net_heat_flux', ...
                'gom_hycom_dt_heat_flux', ...
                'ncep_net_heat_flux', 'monthly_nocs_net_heat_flux', ...
                'landy_net_heat_flux', 'ndbc_sea_t_implied_heat_flux', ...
               };
  legs = {      'TOGA-COARE/NCEP NARR Q_0',    'Smith (1988) Bulk Q_0', ...
                'GoM 4km HYCOM + Q_0', ...
                'NCEP NARR Q_0',      'NOC Southampton Q_0', ...
                'Large & Yeager (2009) Q_0', 'Q_0 "implied" by dT/dt', ...
         };
  for ix = 1:length(fldbasenms)
    flds{ix} = [fldbasenms{ix}];
    mns(ix,1:12) = grpstats(real(stn.(flds{ix}).data),get_month(stn.(flds{ix}).date));
  end;
  figure; maxigraph; hold on;
  plot(1:12,mns); xlim([1 12]);
  % legend(strrep(flds,'_','\_'), 'Location','Best');
  legend(legs, 'Location','Best');
  titlename([stnm ' Net Heat Flux']);
  ylabel('Monthly Mean Flux [ W / m^2 ]');
  % print('-dtiff','../figs/fwyf1-monthly-climatology-comparison-all-fluxes.tiff');


  flds = {}; mns = [];
  fldbasenms = { 'netqf_heat_flux', 'gom_hycom_netqf_heat_flux', ...
                 'ndbc_sea_t_implied_heat_flux' };
  fldbasenms = { 'netqf', 'gom_hycom_netqf', ...
                 'ndbc_sea_t_diff' };
  legs =       { 'HC(Q_0)', 'HC(Q_0 + GoM 4km HYCOM)', ...
                 '\delta_1_h sea temperature' };
  for ix = 1:length(fldbasenms)
    flds{ix} = [fldbasenms{ix}];
    mns(ix,1:12) = grpstats(real(stn.(flds{ix}).data),get_month(stn.(flds{ix}).date));
  end;
  figure; maxigraph; hold on;
  plot(1:12,mns); xlim([1 12]);
  % legend(strrep(flds,'_','\_'), 'Location','Best');
  legend(legs, 'Location','Best');
  titlename([stnm ' Total Heat Budget']);
  ylabel('Mean Temperature Change [ ^oC / hr ]');
  % print('-dtiff','../figs/fwyf1-monthly-climatology-comparison-netqf-fluxes.tiff');

return;
