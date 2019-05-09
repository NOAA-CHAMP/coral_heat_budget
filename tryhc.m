function stn = tryhc(stn,doabs,R,cfac,wfac)
%function stn = tryhc(stn,doabs,R,cfac,wfac)
%
% Do trials of various parameterizations for Horizontal Convection, on the
% heat budget (see STATION_HEAT_BUDGET, CALCQ0, CALCDT) for station STN.  If
% true, optional arg DOABS calculates net absorbed insolation from historical
% KdPAR, site depth (incl. tidal variation), and bottom reflectance for STN.
% Arg R specifies entrainment/mixing (DEFAULT: R=1.00-0.08, 8% mixing eff.)
% CFAC controls proportionality constant for unsteady heat/advective inertia
% scaling (v. Monismith et al. 2006) during cooling episodes (DEFAULT: 1.20);
% WFAC offers the same control for periods of warming (DEFAULT: 0.60).
%
% Last Saved Time-stamp: <Sun 2010-10-03 16:35:22 Eastern Daylight Time gramer>

  figspath = get_thesis_path('../figs');

  if ( ~exist('doabs','var') || isempty(doabs) )
    doabs = false;
  end;
  if ( ~exist('R','var') || isempty(R) )
    % R = 1.00;
    % R = 1.00 - 0.04;
    % R = 1.00 - 0.06;
    R = 1.00 - 0.08;
    % R = 1.00 - 0.11;
    % R = 1.00 - 0.20;
    % R = 1.00 - 0.30;
  end;

  if ( ~exist('cfac','var') || isempty(cfac) )
    cfac = 1.20;
  end;
  if ( ~exist('wfac','var') || isempty(wfac) )
    wfac = 0.60;
  end;

  if ( ~isfield(stn,'tmd_tide') )
    stn = station_tmd_tide(stn,stn.ndbc_sea_t.date);
  end;

  %DEBUG:  hfld = 'tide_i_depth';
  hfld = 'DEBUGtide_i_depth';
  if ( ~isfield(stn,hfld) )
    hfld = 'tmd_tide_i_depth';
    if ( ~isfield(stn,hfld) || all(isnan(stn.(hfld).data)) )
      hfld = 'tpxo_tide_i_depth';
    end;
  end;
  %DEBUG:  hfld = 'tmd_tide_i_depth';
  %DEBUG:  hfld = 'tpxo_tide_i_depth';
  %DEBUG:
  disp(hfld);


  %% TOGA-COARE 3.0a (or 3.0) turbulent, and NCEP NARR radiative, fluxes
  % PFX = 'ndbc_ncep_30_';
  PFX = 'ndbc_ncep_30a_';
  lrfld = 'ncep_lrf';

  % PFX = 'ndbc_ncep_rh_30a_';
  % lrfld = 'ndbc_ncep_rh_lrf';


  %DEBUG:
  lrfld, PFX,


  lffld = [PFX 'latent_heat_flux'];
  sffld = [PFX 'sensible_heat_flux'];
  rffld = [PFX 'rain_heat_flux'];

  if ( ~isfield(stn,lffld) )
    warning('No %s fluxes: trying HFBULK (Smith) instead',PFX);
    %% Simple bulk formulae using Smith (1988) C_d, C_theta, C_q
    PFX = 'ndbc_hfbulk_';
    lffld = [PFX 'latent_heat_flux'];
    sffld = [PFX 'sensible_heat_flux'];
    rffld = '';
  end;
  if ( ~isfield(stn,lffld) )
    error('Which flux?? Neither NDBC_NCEP_* or NDBC_HFBULK found.');
  end;


  % DOABS: Try to estimate total radiation absorbed by bottom/water column?
  % hffld = 'ndbc_ncep_30a_net_heat_flux';
  if (doabs)
    srfld = 'ncep_srf_absorbed';
    hffld = [PFX 'absorbed_heat_flux'];
    ntfld = [PFX 'absorbed_heat_flux_term'];
    dtfld = [PFX 'absorbed_dt'];
%%%% ??? DEBUG
%     if ( ~isfield(stn,srfld) )
      stn = station_absorbed_insolation(stn,srfld,'ncep_srf',hfld);
%     end;
  else
    srfld = 'ncep_srf';
    hffld = [PFX 'total_heat_flux'];
    ntfld = [PFX 'total_heat_flux_term'];
    dtfld = [PFX 'total_dt'];
  end;

  [swix,lwix,lfix,sfix] = intersect_all_dates([], ...
                                              stn.(srfld).date, ...
                                              stn.(lrfld).date, ...
                                              stn.(lffld).date, ...
                                              stn.(sffld).date ...
                                              );
  stn.(hffld).date = stn.(srfld).date(swix);
  stn.(hffld).data = stn.(srfld).data(swix) + ...
      stn.(lrfld).data(lwix) + ...
      stn.(lffld).data(lfix) + ...
      stn.(sffld).data(sfix);

  % Include precipitation flux if we got it
  if ( isfield(stn,rffld) )
    [rfix,hfix] = intersect_dates(stn.(rffld).date,stn.(hffld).date);
    stn.(hffld).date = stn.(hffld).date(hfix);
    stn.(hffld).data = stn.(hffld).data(hfix) + stn.(rffld).data(rfix);
  end;

  % Gross Quality Control - should have been done before now, but...
  badix = find(abs(stn.(hffld).data) > 2500);
  stn.(hffld).date(badix) = [];
  stn.(hffld).data(badix) = [];


  % Calculate heat budget term (Q0/Cp*rho*h)
  stn = station_heat_flux_term(stn,hffld,ntfld,'ndbc_sea_t',[],hfld);


  % Calculate advected heat
  % UdotdelT = 'global_hycom_advected_heat';
  % if ( isfield(stn,'global_hycom_u') )
  %   stn.global_hycom_blank = stn.global_hycom_u;
  %   stn.global_hycom_blank.data(:) = 0;
  %   stn = station_advect_field(stn,UdotdelT,...
  %                              'global_hycom_u','global_hycom_blank',...
  %                              'avhrr_weekly_sst');
  % end;
  UdotdelT = 'stokes_drift_advected_heat';
  if ( isfield(stn,'stokes_drift_u') )
    ucurr = 'stokes_drift_v_7_day_average';
    stn = verify_variable(stn,ucurr);
    stn.stokes_drift_blank = stn.(ucurr);
    stn.stokes_drift_blank.data(:) = 0;
    stn = station_advect_field(stn,UdotdelT,...
                               ucurr,'stokes_drift_blank',...
                               'avhrr_weekly_sst');
  end;
  % % Add heat advection term to total
  % if ( isfield(stn,UdotdelT) )
  %   %DEBUG:
  %   disp([ntfld ' + ' UdotdelT]);
  %   [ix1,ix2] = intersect_dates(stn.(ntfld).date, ...
  %                               stn.(UdotdelT).date);
  %   stn.(dtfld).date = stn.(ntfld).date(ix1);
  %   stn.(dtfld).data = stn.(ntfld).data(ix1) ...
  %       + stn.(UdotdelT).data(ix2);
  % else
  %   %DEBUG:
  %   disp(['No field ' UdotdelT]);
  %   dtfld = hffld;
  % end;


  % Add term for horizontal convection / thermal siphon

  hafld = [hffld '_24_hour_average'];
  % This would cause overdamping of heating and cooling
  % hafld = [hffld];
  % These would eliminate intraday variability from heat flux
  % hafld = [hffld '_24_hour_lowpass'];
  % hafld = [hffld '_40_hour_lowpass'];
  %DEBUG:
  hafld,
  stn = verify_variable(stn,hafld);
  hcfld = ['netqf'];
  % Onset of unstable horizontal convection lags forcing by LAGOFF [hrs]
  lagoff = 0;
  stn = station_horizontal_convection(stn,'ndbc_sea_t',[],hfld,...
                                      hafld,'qvf','qf',R,cfac,wfac,...
                                      ntfld,hcfld,lagoff);

  % Add heat advection term to total
  if ( isfield(stn,UdotdelT) )
    %DEBUG:
    disp([dtfld ' = ' hcfld ' + ' UdotdelT]);
    [ix1,ix2] = intersect_dates(stn.(hcfld).date, ...
                                stn.(UdotdelT).date);
    stn.(dtfld).date = stn.(hcfld).date(ix1);
    stn.(dtfld).data = stn.(hcfld).data(ix1) ...
        + stn.(UdotdelT).data(ix2);
  else
    %DEBUG:
    disp(['No field ' UdotdelT]);
    dtfld = hffld;
  end;

  yrs = [];
  for yr = 1987:2010
    yrix = find(get_year(stn.(dtfld).date) == yr);
    if ( (length(yrix) < 3000) || (max(diff(stn.(dtfld).date(yrix))) > 14) )
      % warning('Too little data to plot %04d',yr);
    else
      yrs(end+1) = yr;
    end;
  end;

  plot_fluxes(stn,yrs(1),1,(yrs(end)-yrs(1)+1)*365);
  datetick('x',2,'keeplimits');
  appendtitlename(sprintf(' (Abs:%g R:%g c:%g w:%g)',doabs,R,cfac,wfac));
%%%%
  % figstub = fullfile(figspath,sprintf('%s_tryhc_%04d-%04d',stn.station_name,yrs(1),yrs(end)));
  % print('-dpng',[figstub '.png']);
%%%%

return;
