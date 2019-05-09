function stn = compare_monthly_flux_climatologies(stn,sfld,fluxfld,doPrint,fnameaddendum,figspath)
%function stn = compare_monthly_flux_climatologies(stn,sfld,fluxfld,doPrint,fnameaddendum,figspath)
%
% Create a plot comparing multiple estimates of monthly mean net heat flux:
% actual (from sea temperature field SFLD), Gramer & Mariano (FLUXFLD), and
% each of ERA-Interim, NCEP NARR, OAFlux/ISCCP, Large&Yeager, NOCS v2. (GoM
% HYCOM NOGAPS forcing is only available over a limited span of 2003-2008.)
%
% If optional DOPRINT, save plot as TIFF in FIGSPATH (DEFAULT: thesis/figs).
%
% Last Saved Time-stamp: <Tue 2013-06-18 14:54:45 Eastern Daylight Time gramer>

  commentstr = '';
  if ( isfield(stn,'commentstr') )
    commentstr = stn.commentstr;
  end;

  if ( ~exist('sfld','var') || isempty(sfld) )
    sfld = 'ndbc_sea_t';
  end;
  if ( ~isfield(stn,sfld) )
    error('No sea temperature field STN.%s',sfld);
  end;

  if ( ~exist('fluxfld','var') || isempty(fluxfld) )
    % fluxfld = 'b_ndbc_erai_erai_30a_net_flux';
    fluxfld = 'simple_ndbc_erai_erai_30a_net_flux';
  end;
  if ( ~isfield(stn,fluxfld) )
    error('No net heat flux field STN.%s',fluxfld);
  end;

  dsfld = [sfld '_diff'];
  if ( ~isfield(stn,dsfld) )
    %stn.(dsfld).date = stn.(sfld).date(1:end-1);
    stn.(dsfld).date = stn.(sfld).date(2:end);
    stn.(dsfld).data = diff(stn.(sfld).data);
    stn = filter_gaps(stn,sfld,dsfld,(1.5/24));
  end;

  dsffld = [dsfld '_flux'];
  if ( ~isfield(stn,dsffld) )
    hfld = ['tpxo_i_depth'];
    mhfld = ['mean_' hfld];
    if ( isfield(stn,mhfld) )
      stn = station_heat_flux_term_inverse(stn,dsffld,dsfld,sfld,[],mhfld);
    else
      warning('No mean water depth field STN.%s found! Using 10m...',mhfld);
      stn = station_heat_flux_term_inverse(stn,dsffld,dsfld,sfld,[],10);
    end;
  end;

  if ( ~exist('doPrint','var') || isempty(doPrint) )
    doPrint = false;
  end;
  if ( ~exist('figspath','var') || isempty(figspath) )
    figspath = get_thesis_path('../figs');
  end;
  if ( ~exist('fnameaddendum','var') || isempty(fnameaddendum) )
    fnameaddendum = '';
  end;

  if ( ~isfield(stn,'erai_net_heat_flux') )
    stn = get_erai_station(stn);
  end;
  if ( ~isfield(stn,'erai_actual_net_heat_flux') )
    stn.erai_turbulent_heat_flux = ts_op(stn.erai_latent_heat_flux,stn.erai_sensible_heat_flux,'+');
    stn.erai_radiative_heat_flux = ts_op(stn.erai_srf,stn.erai_lrf,'+');
    stn.erai_actual_net_heat_flux = ts_op(stn.erai_turbulent_heat_flux,stn.erai_radiative_heat_flux,'+');
  end;
  if ( ~isfield(stn,'ncep_net_heat_flux') )
    stn = get_ncep_station(stn,'narr');
  end;
  if ( ~isfield(stn,'daily_oaflux_net_heat_flux') )
    stn = station_load_oaflux(stn);
  end;
  %if ( ~isfield(stn,'gom_hycom_net_heat_flux') )
  %  stn = get_gom_hycom(stn);
  %end;
  if ( ~isfield(stn,'landy_net_heat_flux') )
    stn = station_load_landy(stn);
  end;
  if ( ~isfield(stn,'monthly_nocs_net_heat_flux') )
    stn = station_load_nocs(stn);
  end;

  if ( ~isfield(stn,'ndbc_hfbulk_net_heat_flux') )
    %stn = station_ndbc_hfbulk(stn,'erai_srf','erai_lrf','erai_relhumid');
    if ( isfield(stn,'erai_ndbc_albedo') )
      alb.date = stn.erai_ndbc_albedo.date;
      alb.data = (1 - stn.erai_ndbc_albedo.data);
    else
      disp([mfilename,' using ERAI albedo for bulk']);
      alb.date = stn.erai_albedo.date;
      alb.data = (1 - stn.erai_albedo.data);
    end;
    stn.bulk_srf = ts_op(alb,stn.erai_dsrf,'.*');
    stn = station_bulk_longwave(stn,'ndbc_air_t','erai_spechumid','erai_barom',...
                                'erai_cloud_cover','ndbc_sea_t','erai_cloud_cover',...
                                'bulk_dlrf','bulk_ulrf','bulk_lrf');
    stn = station_ndbc_hfbulk(stn,'bulk_srf','bulk_lrf','erai_relhumid');
  end;

  % [six,qix,gix,bix,eix,nix,oix,lix,xix] = ...
  %                         stn.gom_hycom_net_heat_flux.date,...
  [six,qix,bix,eix,nix,oix,lix,xix] = ...
      intersect_all_dates([],...
                          stn.(dsffld).date,...
                          stn.(fluxfld).date,...
                          stn.ndbc_hfbulk_net_heat_flux.date,...
                          stn.erai_actual_net_heat_flux.date,...
                          stn.ncep_net_heat_flux.date,...
                          stn.daily_oaflux_net_heat_flux.date,...
                          stn.landy_net_heat_flux.date,...
                          stn.nocs_net_heat_flux.date);

  if ( isempty(six) )
    warning('No dates in "%s" matching all climatologies!',dsffld);
    %%%%%%%%%%%%%%%
    % EARLY RETURN
    %%%%%%%%%%%%%%%
    return;
  end;

  s.date = stn.(dsffld).date(six(1):six(end));                   s.data = stn.(dsffld).data(six(1):six(end));
  q.date = stn.(fluxfld).date(qix(1):qix(end));                  q.data = stn.(fluxfld).data(qix(1):qix(end));
  b.date = stn.ndbc_hfbulk_net_heat_flux.date(bix(1):bix(end));  b.data = stn.ndbc_hfbulk_net_heat_flux.data(bix(1):bix(end));
  %g.date = stn.gom_hycom_net_heat_flux.date(gix(1):gix(end));    g.data = stn.gom_hycom_net_heat_flux.data(gix(1):gix(end));
  e.date = stn.erai_actual_net_heat_flux.date(eix(1):eix(end));  e.data = stn.erai_actual_net_heat_flux.data(eix(1):eix(end));
  n.date = stn.ncep_net_heat_flux.date(nix(1):nix(end));         n.data = stn.ncep_net_heat_flux.data(nix(1):nix(end));
  o.date = stn.daily_oaflux_net_heat_flux.date(oix(1):oix(end)); o.data = stn.daily_oaflux_net_heat_flux.data(oix(1):oix(end));
  l.date = stn.landy_net_heat_flux.date(lix(1):lix(end));        l.data = stn.landy_net_heat_flux.data(lix(1):lix(end));
  x.date = stn.nocs_net_heat_flux.date(xix(1):xix(end));         x.data = stn.nocs_net_heat_flux.data(xix(1):xix(end));
  
  % s.date = stn.(dsffld).date;                   s.data = stn.(dsffld).data;
  % q.date = stn.(fluxfld).date;                  q.data = stn.(fluxfld).data;
  % b.date = stn.ndbc_hfbulk_net_heat_flux.date;  b.data = stn.ndbc_hfbulk_net_heat_flux.data;
  % %g.date = stn.gom_hycom_net_heat_flux.date;    g.data = stn.gom_hycom_net_heat_flux.data;
  % e.date = stn.erai_actual_net_heat_flux.date;  e.data = stn.erai_actual_net_heat_flux.data;
  % n.date = stn.ncep_net_heat_flux.date;         n.data = stn.ncep_net_heat_flux.data;
  % o.date = stn.daily_oaflux_net_heat_flux.date; o.data = stn.daily_oaflux_net_heat_flux.data;
  % l.date = stn.landy_net_heat_flux.date;        l.data = stn.landy_net_heat_flux.data;
  % x.date = stn.nocs_net_heat_flux.date;         x.data = stn.nocs_net_heat_flux.data;

  legs={};
  % legs{end+1}=['Actual (',num2str(get_year(s.date(1))),'-',num2str(get_year(s.date(end))),')'];
  % legs{end+1}=['Gramer&Mariano (',num2str(get_year(q.date(1))),'-',num2str(get_year(q.date(end))),')'];
  % %legs{end+1}=['GoM HYCOM / NOGAPS (',num2str(get_year(g.date(1))),'-',num2str(get_year(g.date(end))),')'];
  % legs{end+1}=['ERA-Interim (',num2str(get_year(e.date(1))),'-',num2str(get_year(e.date(end))),')'];
  % legs{end+1}=['NCEP NARR (',num2str(get_year(n.date(1))),'-',num2str(get_year(n.date(end))),')'];
  % legs{end+1}=['OAFlux / ISCCP 1^o (',num2str(get_year(o.date(1))),'-',num2str(get_year(o.date(end))),')'];
  % legs{end+1}=['Large&Yeager CORE.2 (',num2str(get_year(l.date(1))),'-',num2str(get_year(l.date(end))),')'];
  % legs{end+1}=['NOCS v2 (',num2str(get_year(x.date(1))),'-',num2str(get_year(x.date(end))),')'];
  legs{end+1}=['Actual'];
  legs{end+1}=['Gramer&Mariano Q_0'];
  %legs{end+1}=['Bulk (Smith 1988)'];
  %legs{end+1}=['GoM HYCOM / NOGAPS'];
  legs{end+1}=['ERA-Interim 1.5^o'];
  legs{end+1}=['NCEP NARR 32 km'];
  legs{end+1}=['OAFlux+ISCCP 1^o'];
  legs{end+1}=['L&Y CORE.2'];
  legs{end+1}=['NOCS v2'];

  fmg;
  grpplot_ts(s,@get_month,@nanmean,0,'ks-','LineWidth',2.0);
  grpplot_ts(q,@get_month,@nanmean,0,'ks-.','LineWidth',2.0);
  %grpplot_ts(b,@get_month,@nanmean,0,'k.:');
  %grpplot_ts(g,@get_month,@nanmean,0,'bo:');
  grpplot_ts(e,@get_month,@nanmean,0,'bo-');
  grpplot_ts(n,@get_month,@nanmean,0,'bo-.');
  grpplot_ts(o,@get_month,@nanmean,0,'r^-');
  grpplot_ts(l,@get_month,@nanmean,0,'r^-.');
  grpplot_ts(x,@get_month,@nanmean,0,'r^:');
  %%ylim([-300,150]);
  %ylim([-200,150]);
  ylim([-180,120]);
  lh=legend(legs, 'Location','South');
  %set(lh,'FontSize',14);
  set(lh,'FontSize',15);
  set(gca,'FontSize',14);
  ylabel('W/m^2');
  titlename([upper(stn.station_name) ' Mean Monthly Q_0 ',commentstr,' (',...
             num2str(get_year(s.date(1))),'-',num2str(get_year(s.date(end))),')']);
  set(gca,'color',[.95 .95 .95]);

  if ( doPrint )
    print('-dtiff',fullfile(figspath,[lower(stn.station_name),fnameaddendum,'-',fluxfld,'-',mfilename,'.tiff']));
  end;

  if ( nargout < 1 )
    stn = []; clear stn;
  end;

return;
