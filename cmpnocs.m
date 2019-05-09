function [stn,clim] = cmpnocs(stn,PFX)
%function [stn,clim] = cmpnocs(stn,PFX)
%
% Compare heat budget terms having prefix PFX (DEFAULT: 'ndbc_ncep_30a') with
% monthly National Oceanographic Centre at Southampton v2 heat budget (NOCS,
% 2009; former versions referred to as SOC v1 and v1.1 climatologies).
%
% Last Saved Time-stamp: <Wed 2012-03-28 12:58:48  Lew.Gramer>

  datapath = get_thesis_path('../data');

  if ( ~exist('PFX','var') || ~ischar(PFX) )
    PFX = 'ndbc_ncep_30a';
  end;
  if ( ~isempty(PFX) && PFX(end) ~= '_' )
    PFX(end+1) = '_';
  end;

  if ( ~isfield(stn,'nocs_srf') )
    stn = station_load_nocs(stn);
  end;
  [yr,mo,dy] = datevec(stn.nocs_net_heat_flux.date);
  clim.nocs_srf = grpstats(stn.nocs_srf.data,mo);
  clim.nocs_lrf = grpstats(stn.nocs_lrf.data,mo);
  clim.nocs_latent_heat_flux = grpstats(stn.nocs_latent_heat_flux.data,mo);
  clim.nocs_sensible_heat_flux = grpstats(stn.nocs_sensible_heat_flux.data,mo);
  clim.nocs_net_heat_flux = grpstats(stn.nocs_net_heat_flux.data,mo);
  clim.nocs_heat_flux_term = grpstats(stn.nocs_heat_flux_term.data,mo);

  if ( ~isfield(stn,'ncep_srf') )
    error('No NCEP data in STN!');
    stn = get_ncep_station(stn,'narr');
  end;
  if ( ~isfield(stn,'ncep_latent_heat_flux') )
    x = get_ncep_station(stn.station_name,'narr');
    stn.ncep_latent_heat_flux = x.ncep_latent_heat_flux;
    stn.ncep_sensible_heat_flux = x.ncep_sensible_heat_flux;
    x = []; clear x;
  end;
  [yr,mo,dy] = datevec(stn.ncep_net_heat_flux.date);
  clim.ncep_srf = grpstats(stn.ncep_srf.data,mo);
  clim.ncep_lrf = grpstats(stn.ncep_lrf.data,mo);
  clim.ncep_latent_heat_flux = grpstats(stn.ncep_latent_heat_flux.data,mo);
  clim.ncep_sensible_heat_flux = grpstats(stn.ncep_sensible_heat_flux.data,mo);
  clim.ncep_net_heat_flux = grpstats(stn.ncep_net_heat_flux.data,mo);


  [stn.monthly_sea_t,clim.sea_t,stn.monthly_sea_t_anom] = ...
      monthly_clim_ts(stn.ndbc_sea_t);

  dif.date = stn.monthly_sea_t.date(2:end);
  dif.data = diff(stn.monthly_sea_t.data);

  d = findiff_ts(stn.ndbc_sea_t,7);
  [stn.monthly_sea_t_diff,clim.sea_t_diff,stn.monthly_sea_t_diff_anom] = ...
      monthly_clim_ts(d);

  srfld = 'ncep_srf';
  lrfld = 'ncep_lrf';
  lhfld = [PFX 'latent_heat_flux'];
  shfld = [PFX 'sensible_heat_flux'];
  nhfld = [PFX 'net_heat_flux'];
  ntfld = [PFX 'heat_flux_term'];
  [stn.monthly_srf,clim.srf,stn.monthly_srf_anom] = monthly_clim_ts(stn.(srfld));
  [stn.monthly_lrf,clim.lrf,stn.monthly_lrf_anom] = monthly_clim_ts(stn.(lrfld));
  [stn.monthly_latent_heat_flux,clim.latent_heat_flux,stn.monthly_latent_heat_flux_anom] = monthly_clim_ts(stn.(lhfld));
  [stn.monthly_sensible_heat_flux,clim.sensible_heat_flux,stn.monthly_sensible_heat_flux_anom] = monthly_clim_ts(stn.(shfld));
  [stn.monthly_net_heat_flux,clim.net_heat_flux,stn.monthly_net_heat_flux_anom] = monthly_clim_ts(stn.(nhfld));
  [stn.monthly_heat_flux_term,clim.heat_flux_term,stn.monthly_heat_flux_term_anom] = monthly_clim_ts(stn.(ntfld));
  [stn.monthly_qf,clim.qf,stn.monthly_qf_anom] = monthly_clim_ts(stn.qf);

  % Reproducing graphing done by Qiu et al. (2004)
  clim.nocs_radiative_heat_flux = clim.nocs_srf + clim.nocs_lrf;
  clim.nocs_turbulent_heat_flux = clim.nocs_latent_heat_flux + clim.nocs_sensible_heat_flux;
  clim.ncep_radiative_heat_flux = clim.ncep_srf + clim.ncep_lrf;
  clim.ncep_turbulent_heat_flux = clim.ncep_latent_heat_flux + clim.ncep_sensible_heat_flux;
  clim.radiative_heat_flux = clim.srf + clim.lrf;
  clim.turbulent_heat_flux = clim.latent_heat_flux + clim.sensible_heat_flux;

  % for fldc = {'srf','lrf','latent_heat_flux','sensible_heat_flux','net_heat_flux','qf'}
  for fldc = {'net_heat_flux','qf'}
    fld = fldc{:};
    monfld = ['monthly_' fld];
    nocsfld = ['monthly_nocs_' fld];
    maxigraph(figure);
    hold on;
    plot(stn.(monfld).date,stn.(monfld).data,'-','Color',[.3 .3 .3]);
    plot(stn.(nocsfld).date,stn.(nocsfld).data,'--','Color',[.8 .8 .8]);
    datetick3;
    legend(strrep(PFX,'_','\_'),'NOCS v2');
    title([stn.station_name '.' strrep(fld,'_','\_')]);
  end;


  maxigraph(figure);
  hold on;
  monfld = 'monthly_heat_flux_term';
  nocsfld = 'monthly_nocs_heat_flux_term';
  if ( ~isfield(stn,nocsfld) )
    stn = station_heat_flux_term(stn,'monthly_nocs_net_heat_flux',...
                                 nocsfld,'monthly_nocs_sea_t',[],[]);
  end;

  plot(stn.(monfld).date,stn.(monfld).data,'-','Color',[.3 .3 .3]);
  plot(stn.(nocsfld).date,stn.(nocsfld).data,'--','Color',[.7 .7 .7]);
  plot(dif.date,dif.data,':','Color',[.0 .0 .0]);
  datetick3;
  legend(strrep(PFX,'_','\_'),'NOCS v2','\Delta(T_s)/\Delta t');
  title([stn.station_name '.' strrep('heat_flux_term','_','\_')]);


  % maxigraph(figure);
  % hold on;
  % monfld = 'monthly_heat_flux_term';
  % nocsfld = 'monthly_nocs_heat_flux_term';

  % plot(stn.(monfld).date,stn.(monfld).data,'-','Color',[.3 .3 .3]);
  % plot(stn.(nocsfld).date,stn.(nocsfld).data,'--','Color',[.7 .7 .7]);
  % plot(stn.monthly_sea_t_diff.date,stn.monthly_sea_t_diff.data,':','Color',[.0 .0 .0]);
  % datetick3;
  % legend(strrep(PFX,'_','\_'),'NOCS v2','\Delta(T_s)/\Delta(t)');
  % title([stn.station_name '.' strrep('heat_flux_term','_','\_')]);


  % figure; subplot(2,1,1);plot(1:12,clim.nocs_turbulent_heat_flux,'o-',1:12,clim.ncep_turbulent_heat_flux,'s-',1:12,clim.turbulent_heat_flux,'.-');xlim([0.5,12.5]);set(gca,'xtick',[1:12],'xticklabel',{'J','F','M','A','M','J','J','A','S','O','N','D'});legend('NOCS','NCEP','3.0a');title('Q_L_H+Q_S_H'); subplot(2,1,2);plot(1:12,clim.nocs_srf+clim.nocs_lrf,'o-',1:12,clim.ncep_srf+clim.ncep_lrf,'s-',1:12,clim.ncep_srf+clim.ncep_lrf,'.-');xlim([0.5,12.5]);set(gca,'xtick',[1:12],'xticklabel',{'J','F','M','A','M','J','J','A','S','O','N','D'});legend('NOCS','NCEP','Bulk');title('Q_S_W+Q_L_W'); maxigraph;

  figure; subplot(2,1,1);plot(1:12,clim.nocs_latent_heat_flux,'o-',1:12,clim.ncep_latent_heat_flux,'s-',1:12,clim.latent_heat_flux,'.-');xlim([0.5,12.5]);set(gca,'xtick',[1:12],'xticklabel',{'J','F','M','A','M','J','J','A','S','O','N','D'});legend('NOCS','NCEP','3.0a');title('Q_L_H'); subplot(2,1,2);plot(1:12,clim.nocs_sensible_heat_flux,'o-',1:12,clim.ncep_sensible_heat_flux,'s-',1:12,clim.sensible_heat_flux,'.-');xlim([0.5,12.5]);set(gca,'xtick',[1:12],'xticklabel',{'J','F','M','A','M','J','J','A','S','O','N','D'});legend('NOCS','NCEP','3.0a');title('Q_S_H'); maxigraph;

  figure; plot(1:12,clim.nocs_net_heat_flux,'o-',1:12,clim.ncep_net_heat_flux,'s-',1:12,clim.net_heat_flux,'.-',1:12,clim.qf,'^-');xlim([0.5,12.5]);set(gca,'xtick',[1:12],'xticklabel',{'J','F','M','A','M','J','J','A','S','O','N','D'});legend('NOCS','NCEP','3.0a','3.0a+HC');title('Q_0'); maxigraph;

return;
