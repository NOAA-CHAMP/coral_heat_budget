function cs = cmpclims(stn,dltflds,accflds,climfld,meanfld)
%function cs = cmpclims(stn,dltflds,accflds,climfld,meanfld)

  if ( ~exist('dltflds','var') || isempty(dltflds) )
    %dltflds={'ndbc_sea_t','ndbc_air_sea_t','ndbc_air_sea_spechumid'}; accflds={'qswterm','ndbc_bulk_rh_lrf_term','ndbc_bulk_30_sensible_flux_term','ndbc_bulk_30_latent_flux_term','q0term'}; ttl='NOCS, NCEP, Bulk, TOGA-COARE 2.0, 2.6, 3.0';

    %dltflds={'ndbc_sea_t','ndbc_air_sea_t','ndbc_air_sea_spechumid'}; accflds={'nocs_heat_flux_term','ncep_heat_flux_term','ncep_bulk_heat_flux_term','ndbc_hfbulk_heat_flux_term','ndbc_bulk_30_heat_flux_term'}; ttl='NOCS, NCEP, NCEP/Bulk, Bulk, TOGA-COARE 2.0, 2.6, 3.0';

    %dltflds={'ndbc_sea_t'}; accflds={'nocs_heat_flux_term','ncep_heat_flux_term','ncep_20_heat_flux_term','ndbc_hfbulk_heat_flux_term','ndbc_bulk_rh_heat_flux_term','ndbc_bulk_26_heat_flux_term','ndbc_bulk_30_heat_flux_term','ndbc_bulk_30a_heat_flux_term'}; ttl='NOCS, NCEP, NCEP/Bulk, Bulk, TOGA-COARE 2.0,2.6,3.0,3.0a';
    %dltflds={'ndbc_sea_t'}; accflds={'ncep_heat_flux_term','ncep_bulk_heat_flux_term','ndbc_hfbulk_heat_flux_term','ndbc_bulk_rh_heat_flux_term','ndbc_bulk_26_heat_flux_term','ndbc_bulk_30_heat_flux_term','ndbc_bulk_30a_heat_flux_term','ndbc_ncep_30a_heat_flux_term'}; ttl='NCEP, NCEP/Bulk, Bulk, TOGA-COARE 2.6,3.0,3.0a,3.0a+NCEP LW';
    %dltflds={'ndbc_sea_t'}; accflds={'ncep_heat_flux_term','ndbc_ncep_26_heat_flux_term','ndbc_hfbulk_heat_flux_term','ndbc_bulk_rh_heat_flux_term','ndbc_bulk_26_heat_flux_term','ndbc_bulk_30_heat_flux_term','ndbc_bulk_30a_heat_flux_term','ndbc_ncep_30a_heat_flux_term'}; ttl='NCEP, NCEP/Bulk, Bulk, TOGA-COARE 2.6,3.0,3.0a,3.0a+NCEP LW';

    %dltflds={'ndbc_sea_t','ndbc_wind1_speed'}; accflds={'ncep_heat_flux_term','ndbc_ncep_26_heat_flux_term','ndbc_ncep_30_heat_flux_term','ndbc_ncep_30a_heat_flux_term'}; ttl='NCEP vs. NCEP/TOGA-COARE 2.6,3.0,3.0a';

    %dltflds={'ndbc_sea_t','ndbc_wind1_speed'}; accflds={'ncep_heat_flux_term','ndbc_ncep_26_heat_flux_term','ndbc_ncep_30_heat_flux_term','ndbc_ncep_30a_heat_flux_term','ndbc_ncep_30a_dt'}; ttl='NCEP vs. NCEP/TOGA-COARE 2.6,3.0,3.0a,dT/dt';
    dltflds={'ndbc_sea_t'};
    accflds={'nocs_heat_flux_term','ncep_heat_flux_term','ndbc_ncep_30a_heat_flux_term','ndbc_ncep_30a_dt','netqf'};
    ttl='NOCS, NCEP, NCEP/TOGA-COARE 3.0a with and w/o advection, and dT/dt';

    %dltflds={'ndbc_sea_t','ndbc_air_sea_t','ndbc_air_sea_spechumid'}; accflds={'nocs_heat_flux_term','ncep_heat_flux_term','ndbc_rhbulk_heat_flux_term','ndbc_bulk_rh_heat_flux_term','ndbc_bulk_26_heat_flux_term','ndbc_bulk_30_heat_flux_term'}; ttl='NOCS, NCEP, Bulk, TOGA-COARE 2.0, 2.6, 3.0';

    %dltflds={'ndbc_sea_t','ndbc_air_sea_spechumid'}; accflds={'nocs_latent_flux_term','ncep_latent_flux_term','ndbc_rhbulk_latent_flux_term','ndbc_bulk_rh_latent_flux_term','ndbc_bulk_26_latent_flux_term','ndbc_bulk_30_latent_flux_term'}; ttl='Q_L_H: NOCS, NCEP, Bulk, TOGA-COARE 2.0, 2.6, 3.0';

    %dltflds={'ndbc_sea_t'}; accflds={'ncep_heat_flux_term','ncep_lrf_term','ncep_latent_flux_term','ncep_sensible_flux_term'}; ttl='ncep';

    %dltflds={'ndbc_sea_t'}; accflds={'ncep_bulk_heat_flux_term','ncep_lrf_term','ncep_bulk_latent_flux_term','ncep_bulk_sensible_flux_term'}; ttl='ncep_bulk';
    %dltflds={'ndbc_sea_t'}; accflds={'ncep_bulk_rh_heat_flux_term','ndbc_bulk_rh_lrf_term','ncep_bulk_rh_latent_flux_term','ncep_bulk_rh_sensible_flux_term'}; ttl='NCEP 2.0';
    %dltflds={'ndbc_sea_t'}; accflds={'ndbc_bulk_heat_flux_term','ndbc_bulk_lrf_term','ndbc_bulk_latent_flux_term','ndbc_bulk_sensible_flux_term'}; ttl='2.0';

    %dltflds={'ndbc_sea_t'}; accflds={'ndbc_bulk_rh_heat_flux_term','ndbc_bulk_rh_lrf_term','ndbc_bulk_rh_latent_flux_term','ndbc_bulk_rh_sensible_flux_term'}; ttl='2.0';
    %dltflds={'ndbc_sea_t'}; accflds={'ndbc_bulk_26_heat_flux_term','ndbc_bulk_rh_lrf_term','ndbc_bulk_26_latent_flux_term','ndbc_bulk_26_sensible_flux_term'}; ttl='2.6';
    %dltflds={'ndbc_sea_t'}; accflds={'ndbc_bulk_30_heat_flux_term','ndbc_bulk_rh_lrf_term','ndbc_bulk_30_latent_flux_term','ndbc_bulk_30_sensible_flux_term','ndbc_bulk_30_rain_flux_term'}; ttl='3.0';
  end;

  if ( ~exist('climfld','var') )
    climfld = 'hrclim';
    %climfld = 'dyclim';
  end;
  if ( ~exist('meanfld','var') )
    meanfld = 'datmtx';
    %meanfld = 'dymean';
  end;

  fldnms(1:length(dltflds)) = strrep(dltflds,'_','\_');
  fldnms(end+1:end+length(accflds)) = strrep(accflds,'_','\_');

  % linspec = {'b.', 'g.', 'r.', 'c.', 'k.', 'y.', 'm.', 'bd', 'gd', 'rd', 'cd', 'kd', 'yd', 'md'};
  % linspec = {'k.-','k.:','k.--','ko-','ko:','ko--','kp-','kp:','kp--','k^-','k^:','k^--','k*-','k*:','k*--'};
  linspec = {'k.-','k.:','ko-','ko:','kp-','kp:','k^-','k^:','k*-','k*:','kd-','kd:','ks-','ks:'};

  cix = 1;
  for ix = 1:length(dltflds)
    cs(cix) = anclim(stn,dltflds{ix});
    cix = cix + 1;
  end;
  for ix = 1:length(accflds)
    if ( ~isfield(stn,accflds{ix}) )
      stn.(accflds{ix}).date = stn.(dltflds{1}).date;
      stn.(accflds{ix}).data = repmat(0,size(stn.(dltflds{1}).data));
    end;
    cs(cix) = anclim(stn,accflds{ix});

    peryr = length(cs(cix).(climfld));
    switch (peryr)
     case 12,   yearly = 1:12;
     case 52,   yearly = 1:52;
     case 365,  yearly = 1:365;
     otherwise, yearly = 0:(1/24):365-(1/48);
    end;
    cix = cix + 1;
  end;

  cix = 1;
  for ix = 1:length(dltflds)
    cs(cix).name = dltflds{ix};
    cix = cix + 1;
  end;
  for ix = 1:length(accflds)
    cs(cix).name = accflds{ix};
    cix = cix + 1;
  end;

  % First compare climatologies
  cix = 1;
  x = yearly;
  for ix = 1:length(dltflds)
    y(cix,1:length(cs(cix).(climfld))) = cs(cix).(climfld) - cs(cix).(climfld)(1);
    cix = cix + 1;
  end;
  for ix = 1:length(accflds)
    y(cix,:) = cumsum(cs(cix).(climfld));
    cix = cix + 1;
  end;

  if (1)
    figure;
    % plot(x,y);
    hold on;
    for cix = 1:size(y,1)
      plot(x,y(cix,:),linspec{mod(cix-1,length(linspec))+1});
    end;
    maxigraph;
    legend(fldnms,'Location','NW');
    titlename(sprintf('Climatologies (%s) - %s',climfld,ttl));
    hold off;
    ylim([-7 14]);
  end;

  figure;
  hold on;
  for cix = 1:length(cs)
    x = [];
    y = [];
    for yix = 1:length(cs(cix).yrs)
      newyear = datenum(cs(cix).yrs(yix),1,1) - 1;
      dts = (yearly + newyear);
      rawdat = cs(cix).(meanfld)(yix,:);
      gdix = find(isfinite(rawdat));
      if ( ~isempty(gdix) )
        dat = interp1(dts(gdix),rawdat(gdix),dts,'linear',0);
      else
        dat = rawdat;
      end;

      x = [ x dts ];
      if ( cix <= length(dltflds) )
        if (dat(1)==0); dat(1) = min(dat(gdix)); end;
        y = [ y dat-dat(1) ];
      else
        y = [ y cumsum(dat) ];
      end;
    end;
    plot(x,y,linspec{mod(cix-1,length(linspec))+1});
  end;
  maxigraph;
  %%%% ??? HACK
  % xlim([datenum(1998,1,1) datenum(2002,2,1)]);
  % xlim([datenum(1996,1,1) datenum(2004,1,1)]);
  ylim([-25 15]);
  datetick3;
  grid on;
  legend(fldnms,'Location','SW');
  titlename(sprintf('Time series (%s) - %s',meanfld,ttl));

  if ( nargout < 1 )
    cs = [];
  end;

return;
