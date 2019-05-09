1;

%perfun = @get_hour;
%perfun = @get_jday;
perfun = @get_yearday;

%stnms = {'fwyf1','mlrf1','lonf1','smkf1','sanf1','dryf1'};
stnms = {'fwyf1','mlrf1','lonf1','smkf1','sanf1'};

%'wind_speed' 'U' [3,17]
vars = {'air_t',    'wind_u',    'wind_v',    'barom'};
ylbs = {'T_a',      'U^x',       'U^y',       'p_a'};
ylms = {[18,32],    [-12,12],    [-12,12],    [1012,1023]};

%avgper = '_6_h_lp';
avgper = '';

if ( ~exist('stns','var') )
    for stix=1:numel(stnms)
        stns{stix} = get_station_from_station_name(stnms{stix});
        stns{stix} = load_all_ndbc_data(stns{stix});
        stns{stix} = get_erai_station(stns{stix});
        stns{stix} = get_ncep_station(stns{stix},'narr');
    end;
end;

disp('Blue:  ERA-Interim');
disp('Red:   NCEP NARR');
disp('Black: NDBC in situ');


nstns = numel(stns);
nvars = numel(vars);

fmg;
for stix=1:nstns
    disp(stns{stix}.station_name);

    for varix=1:nvars
        var = vars{varix};
        ivar = ['ndbc_' var avgper];
        evar = ['erai_' var '_36_h_lp'];
        nvar = ['ncep_' var avgper];
        if ( ~isempty(strfind(ivar,'wind_')) )
            ivar = strrep(ivar,'wind_','wind1_');
        end;
        stns{stix} = verify_variable(stns{stix},ivar);
        stns{stix} = verify_variable(stns{stix},evar);
        stns{stix} = verify_variable(stns{stix},nvar);

        [iix,eix,nix] = intersect_all_dates([],...
                                            stns{stix}.(ivar).date,...
                                            stns{stix}.(evar).date,...
                                            stns{stix}.(nvar).date);
        [icum,itid]=grp_ts(stns{stix}.(ivar).data(iix),stns{stix}.(ivar).date(iix),perfun);
        [ecum,etid]=grp_ts(stns{stix}.(evar).data(eix),stns{stix}.(evar).date(eix),perfun);
        [ncum,ntid]=grp_ts(stns{stix}.(nvar).data(nix),stns{stix}.(nvar).date(nix),perfun);

        plotix = (varix-1)*nstns+stix;
        subplot_tight(nvars,nstns,plotix);
        plot(etid,ecum,'b',ntid,ncum,'r',itid,icum,'k');
        xlim([min(itid(:))-1,max(itid(:))+1]);
        ylim(ylms{varix});
        if (plotix<=nstns); title(upper(stns{stix}.station_name)); end;
        if (mod(plotix,nstns)==1); ylabel(ylbs{varix}); end;

        clear var ivar evar nvar icum itid ecum etid ncum ntid plotix
    end;
end;

if ( 0 )
  % for ix=1:nstns;
  ix=2;
    fh=fmg;
    subplot_tight(2,2,1); scatter_fit_ts(stns{ix}.erai_air_t,stns{ix}.ndbc_air_t,[],@ts_boreal_cool,'ERAI T_a','NDBC T_a',fh,[],true); axis([3,33,3,33]); titlename(stns{ix}.station_name);
    subplot_tight(2,2,2); scatter_fit_ts(stns{ix}.ncep_air_t,stns{ix}.ndbc_air_t,[],@ts_boreal_cool,'NARR T_a','NDBC T_a',fh,[],true); axis([3,33,3,33]); titlename(stns{ix}.station_name);
    subplot_tight(2,2,3); scatter_fit_ts(stns{ix}.erai_wind_speed,stns{ix}.ndbc_wind1_speed,[],@ts_boreal_cool,'ERAI U','NDBC U',fh,[],true); axis([0,60,0,60]); titlename(stns{ix}.station_name);
    subplot_tight(2,2,4); scatter_fit_ts(stns{ix}.ncep_wind_speed,stns{ix}.ndbc_wind1_speed,[],@ts_boreal_cool,'NARR U','NDBC U',fh,[],true); axis([0,60,0,60]); titlename(stns{ix}.station_name);
  % end;
end;
