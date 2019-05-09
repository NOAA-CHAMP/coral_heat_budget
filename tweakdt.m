function stn = tweakdt(stn)
%function stn = tweakdt(stn)

  hfld = 'tmd_tide_i_depth';

%       'ndbc_bulk_rh','ndbc_bulk_26','ndbc_bulk_30','ndbc_bulk_30a',...

  for q0basename = { ...
      'ndbc_hfbulk',...
      'ndbc_ncep_26','ndbc_ncep_30','ndbc_ncep_30a',...
      'bic_ncep','bic_ncep_rh','bic_ncep_26','bic_ncep_30','bic_ncep_30a',...
                   }

    nffld = [q0basename{:} '_net_heat_flux'];
    htfld = [q0basename{:} '_heat_flux_term'];

    if ( ~isfield(station,nffld) )
        %DEBUG:
        disp(['No field ' nffld]);

    else

      if ( ~isempty(strfind('bic',nffld)) )
        srfld = 'bic_surf_srf';
      else
        srfld = 'ncep_srf';
      end;
      safld = [srfld '_absorbed'];

      stn = station_absorbed_insolation(stn,safld,srfld,hfld);

      lrfld = 'ncep_lrf';
      lffld = [q0basename{:} '_latent_heat_flux'];
      sffld = [q0basename{:} '_sensible_heat_flux'];
      udfld = 'advected_heat';

      [saix,lrix,lfix,sfix] = ...
          intersect_all_dates(stn.(safld).date, stn.(lrfld).date, ...
                              stn.(lffld).date, stn.(sffld).date);

      stn.(nffld).date = stn.(safld).date(saix);
      stn.(nffld).data = stn.(safld).data(saix) ...
          + stn.(lrfld).data(lrix) ...
          + stn.(lffld).data(lfix) ...
          + stn.(sffld).data(sfix);

      stn = station_heat_flux_term(stn,nffld,htfld,'ndbc_sea_t',[],hfld);

      qvfld = [q0basename{:} '_thermal_discharge'];
      for R = 0.1:0.05:0.3
        mixstr = sprintf('%03.0f',100*R);
        dtfld = [q0basename{:} '_dt_' mixstr];
        hcfld = [q0basename{:} '_horzconv_' mixstr];
        stn = station_horizontal_convection(stn,'ndbc_sea_t',[],hfld,nffld,qvfld,hcfld,R)

        [htix,udix,hcix] = ...
            intersect_all_dates(stn.(htfld).date,stn.(udfld).date,stn.(hcfld).date);

        stn.(dtfld).date = stn.(htfld).date;
        stn.(dtfld).data = stn.(htfld).data ...
            + stn.(udfld).data(udix) ...
            + stn.(hcfld).data(hcix);
      end;

    end;

  end;

return;
