function stn = anfluxbox(stn)

  % Interannual, monthly, weekly, diurnal
  doPlot = [1,1,1,0];
  doPrint = [1,1,1,0];
  % doPrint = [0,0,0,0];


  % station_boxplots(stn,'ndbc_sea_t_erai_erai_30a_latent_flux','Q_L_H',[-600,100],[],[],0,doPlot,doPrint);
  % station_boxplots(stn,'ndbc_sea_t_erai_erai_30a_sensible_flux','Q_S_H',[-600,100],[],[],0,doPlot,doPrint);
  % station_boxplots(stn,'ndbc_sea_t_erai_erai_30a_turbulent_flux','Q_L_H+Q_S_H+Q_R_H',[-800,100],[],[],0,doPlot,doPrint);
  % station_boxplots(stn,'erai_srf','Q_S_W',[0,600],[],[],0,doPlot,doPrint);
  % station_boxplots(stn,'ndbc_sea_t_erai_erai_lrf','Q_L_W',[-600,100],[],[],0,doPlot,doPrint);
  % station_boxplots(stn,'ndbc_sea_t_erai_erai_rf','Q_S_W+Q_L_W',[-600,1000],[],[],0,doPlot,doPrint);
  % station_boxplots(stn,'ndbc_sea_t_erai_erai_arf','\gammaQ_S_W+Q_L_W',[-600,1000],[],[],0,doPlot,doPrint);
  % station_boxplots(stn,'ndbc_sea_t_erai_erai_30a_net_flux','Q_0',[-1000,1000],[],[],0,doPlot,doPrint);

  % ylm = [-1000,1000];
  ylm = [-1500,1500];

  station_boxplots(stn,'ndbc_sea_t_erai_erai_30a_latent_flux','Q_L_H',ylm,[],[],0,doPlot,doPrint);
  station_boxplots(stn,'ndbc_sea_t_erai_erai_30a_sensible_flux','Q_S_H',ylm,[],[],0,doPlot,doPrint);
  station_boxplots(stn,'ndbc_sea_t_erai_erai_30a_turbulent_flux','Q_L_H+Q_S_H+Q_R_H',ylm,[],[],0,doPlot,doPrint);
  station_boxplots(stn,'erai_srf','Q_S_W',ylm,[],[],0,doPlot,doPrint);
  station_boxplots(stn,'ndbc_sea_t_erai_erai_lrf','Q_L_W',ylm,[],[],0,doPlot,doPrint);
  station_boxplots(stn,'ndbc_sea_t_erai_erai_rf','Q_S_W+Q_L_W',ylm,[],[],0,doPlot,doPrint);
  station_boxplots(stn,'ndbc_sea_t_erai_erai_arf','\gammaQ_S_W+Q_L_W',ylm,[],[],0,doPlot,doPrint);
  station_boxplots(stn,'ndbc_sea_t_erai_erai_30a_net_flux','Q_0',ylm,[],[],0,doPlot,doPrint);
  station_boxplots(stn,'b_erai_qbo','Q_b',ylm,[],[],0,doPlot,doPrint);
  station_boxplots(stn,'b_ndbc_sea_t_erai_erai_30a_net_flux','Q_0+Q_b',ylm,[],[],0,doPlot,doPrint);
  station_boxplots(stn,'ndbc_sea_t_erai_erai_30a_erai_avhrr_hc_dTdthc_flux','\rhoC_ph(u_h_c^.\nabla_hT_h_c)',ylm,[],[],0,doPlot,doPrint);
  station_boxplots(stn,'ndbc_sea_t_erai_erai_30a_erai_avhrr_hc_dTdt_flux','\rhoC_ph\partial_tT_s',ylm,[],[],0,doPlot,doPrint);

  if ( nargout < 1 )
    stn=[]; clear stn;
  end;

return;
