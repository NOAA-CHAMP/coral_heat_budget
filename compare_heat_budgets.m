function stn = compare_heat_budgets(stn,sfld,RAPFX,KMPFX)
%function stn = compare_heat_budgets(stn,sfld,RAPFX,KMPFX)

  if ( ~exist('sfld','var') || isempty(sfld) )
    sfld = 'ndbc_sea_t';
  end;

  if ( ~exist('RAPFX','var') || isempty(RAPFX) )
    RAPFX = 'erai';
    % RAPFX = 'ncep';
  end;

  if ( ~exist('KMPFX','var') || isempty(KMPFX) )
    KMPFX = 'fkeys_hycom';
    % KMPFX = 'gom_hycom';
  end;

  TURPFX = ['ndbc_' RAPFX '_30a'];

  WAVEPFX = 'ww3';

  switch (KMPFX),
   case 'fkeys_hycom',	QEPFX = [WAVEPFX '_fkeys_qe'];
   case 'gom_hycom',	QEPFX = [WAVEPFX '_gom_qe'];
   otherwise,		error('Unknown km-scale model "%s"',KMPFX);
  end;

  % Net surface heat flux
  qtfld = [TURPFX '_net_heat_flux_term'];

  % Total budget
  dTfld = [TURPFX '_' QEPFX '_dt'];

  % Water-benthos fluxes
  qbfld = ['benthic_' RAPFX '_srf'];
  btfld = ['benthic_' RAPFX '_t'];
  qbofld = ['benthic_' RAPFX '_qbo'];

  bdTfld = ['benthic_' dTfld];
  bdTffld = [bdTfld '_heat_flux'];

  [firstyr,ig,ig] = datevec(stn.(dTfld).date(1));
  [lastyr,ig,ig] = datevec(stn.(dTfld).date(end));
  dys = datenum(lastyr,12,31) - datenum(firstyr,1,1) + 1;

  plot_fluxes(stn,firstyr,1,dys,{sfld,btfld},[],{'ndbc_hfbulk_heat_flux_term',qtfld,dTfld,bdTfld,[bdTfld '_netqf']},[],...
              {'NDBC sea temperature','Modeled substrate temperature',...
               'Bulk Q_0',...
               '(Q_0 == \gammaQ_S_W + Q_L_W + Q_L_H + Q_S_H)/\rhoC_ph',...
               'K_\theta_H\nabla^2T_1_k_m + (u_1_k_m + u_q_e)^.\nablaT_1_k_m + Q_0/\rhoC_ph',...
               'K_\theta_H\nabla^2T_1_k_m + (u_1_k_m + u_q_e)^.\nablaT_1_k_m + (Q_0+Q_b)/\rhoC_ph',...
               'HC( K_\theta_H\nabla^2T_1_k_m + (u_1_k_m + u_q_e)^.\nablaT_1_k_m + (Q_0+Q_b)/\rhoC_ph )',...
              });
  appendtitlename(strrep(sprintf(' (%s,%s)', upper(RAPFX), upper(KMPFX)),'_','\_'));

  if ( nargout < 1 )
    stn = []; clear stn;
  end;

return;
