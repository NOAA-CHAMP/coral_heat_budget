function stns = xanndbc(stns)

  datapath = get_thesis_path('../data');
  figspath = get_thesis_path('../figs');
  taufld = 'ndbc_sea_t_erai_erai_30a_wind_stress';

  yrs=1996:2010;

  stnms = { 'fwyf1','mlrf1','lonf1','smkf1','sanf1','dryf1' };

  vars = { 'ndbc_air_t', 'ndbc_sea_t', 'ndbc_wind1_speed', 'tau_xshore', };
  ylbs = { 'T_a',        'T_s',        'U',                '\tau^x^s',   };
  ylms = { [ 4,35],      [ 4,35],      [0,45],             [-0.11,0.11],   };
  % vars = { 'ndbc_air_t', 'ndbc_sea_t', 'ndbc_wind1_speed', 'ndbc_wind1_u', };
  % ylbs = { 'T_a',        'T_s',        'U',                '\U^x^s',   };
  % ylms = { [ 4,34],       [ 4,34],       [0,45],           [-30,+30],   };

  if ( ~exist('stns','var') )
    stns = {};
  end;
  for stix=1:numel(stnms)
    stnm = lower(stnms{stix});

    if ( numel(stns) < stix || ...
         ~isfield(stns{stix},'station_name') || ...
         ~strcmpi(stnms{stix},stns{stix}.station_name) )

      stns{stix} = get_station_from_station_name(stnms{stix});
      stns{stix} = load_all_ndbc_data(stns{stix});

      stns{stix} = verify_variable(stns{stix},'ndbc_wind1_u');
      stns{stix} = verify_variable(stns{stix},'ndbc_wind1_v');
      stns{stix}.ndbc_wind1_u.data = kts2mps(stns{stix}.ndbc_wind1_u.data);
      stns{stix}.ndbc_wind1_v.data = kts2mps(stns{stix}.ndbc_wind1_v.data);

      % Get a very limited subset of the (Mbs worth of) ERA-I reanalysis data
      x = get_erai_station(stnm);
      stns{stix}.erai_air_t = x.erai_air_t;
      stns{stix}.erai_barom = x.erai_barom;
      stns{stix}.erai_relhumid = x.erai_relhumid;
      stns{stix}.erai_wind_speed = x.erai_wind_speed;
      stns{stix}.erai_wind_stress_u = x.erai_wind_stress_u;
      stns{stix}.erai_wind_stress_v = x.erai_wind_stress_v;
      x=[]; clear x;

      stns{stix}.erai_wind_stress.date = stns{stix}.erai_wind_stress_u.date;
      stns{stix}.erai_wind_stress.data = ...
          uv_to_spd(stns{stix}.erai_wind_stress_u.data,...
                    stns{stix}.erai_wind_stress_v.data);

      matfname = fullfile(datapath,[stnm '_ndbc_sea_t_erai_erai_30a_wind_stress.mat']);
      if ( exist(matfname,'file') )
        stns{stix}.tau = load(matfname);
        [tix,dix] = intersect_dates(stns{stix}.tau.date,stns{stix}.ndbc_wind1_dir.date);

        t.date = stns{stix}.tau.date(tix); t.data = stns{stix}.tau.data(tix);
        d.date = stns{stix}.ndbc_wind1_dir.date(dix); d.data = stns{stix}.ndbc_wind1_dir.data(dix);
        stns{stix}.tau_x.date = t.date;
        stns{stix}.tau_y.date = t.date;
        [stns{stix}.tau_x.data,stns{stix}.tau_y.data] = spddir_to_uv(t.data,d.data);

      else
        disp([stnm ' Bulk Windstress']);
        stns{stix} = ...
            station_bulk_windstress(stns{stix},'tau',...
                                    'ndbc_wind1_speed',[],'ndbc_wind1_dir',...
                                    'ndbc_air_t','erai_relhumid','erai_barom');
      end;

      stns{stix} = station_optimal_isobath_orientation(stns{stix});
      stns{stix} = station_reorient_vectors(stns{stix},'isobath_orientation',...
                                            'ndbc_wind1_u','ndbc_wind1_v');
      stns{stix} = station_reorient_vectors(stns{stix},'isobath_orientation',...
                                            'tau_x','tau_y');
      stns{stix} = station_reorient_vectors(stns{stix},'isobath_orientation',...
                                            'erai_wind_stress_u','erai_wind_stress_v');

    end;
  end;


  nstns = numel(stns);
  nvars = numel(vars);

  for varix=1:nvars
    for stix=1:nstns
      % Default: Use ALL available data for each variable...
      ixen{stix,varix} = 1:numel(stns{stix}.(vars{varix}).date);
      % % Use only data from intermediate years
      % ixen{stix,varix} = find(ismember(get_year(stns{stix}.(vars{varix}).date),yrs));
      % % Use only dates when all stations agree!
      % dts{stix} = stns{stix}.(vars{varix}).date;
    end;
    % % Use only dates when all stations agree!
    % res = intersect_all_dates([],dts{:});
    % for stix=1:nstns
    %   ixen{stix,varix} = res{stix};
    % end;
  end;


  ncols = ceil(nstns/2);
  nrows = ceil(nvars/2)*2;

  panelix = 0;
  for vargpix=1:2:nvars
    fmg;
    for varsubix=1:2
      varix = (vargpix-1)+varsubix;
      for stgpix=1:3:nstns
        for stsubix=1:3
          stix=(stgpix-1)+stsubix;
          %plotix = (stix-1)*nvars+varix;
          plotix = (varsubix-1)*(ncols*2)+stix;
          ax=subplot_tight(nrows,ncols,plotix);
          mos = get_month(stns{stix}.(vars{varix}).date(ixen{stix,varix}));
          dat = stns{stix}.(vars{varix}).data(ixen{stix,varix});
          boxplot(dat, mos, 'notch','on', 'whisker',2);
          xlim([0.5 12.5]);
          ylim(ylms{varix});
          xlabel([]); ylabel([]);
          grid on;
          set(ax,'fontsize',10, 'units','normalized');
          set(ax,'xtick',2:2:12,'xticklabel',num2str([2:2:12]'));
          % th = title( ['(' underline(char('a'+(plotix-1)),12) ')'] );
          th = title( ['(' char('a'+panelix+(plotix-1)) ')'] );
          pos = get(th,'position');
          pos(1) = -0.50;
          set(th,'position',pos,'verticalalignment','bottom','fontsize',12,'fontweight','demi');
          if (0)
            % if ( plotix<=nstns ); titlename(upper(stnms{stix})); end;
            title(upper(stnms{stix}));
            if ( mod(plotix,nstns) == 1 ); ylabel(ylbs{varix}); end;
            % ylabel(ylbs{varix}); 
            %xlabel('Year-Month');
            %ylabel([upper(stnms{stix}) ' ' ylbs{varix}]);
          end;
        end;
      end;
    end;
    print('-dtiff',fullfile(figspath,['fig3-panel-' num2str(panelix) '.tiff']));
    panelix = panelix + (nrows*ncols);
  end;

return;
