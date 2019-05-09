1;

diary(fullfile(get_thesis_path('../src'),'reopt_all.log'));
more off;

for cstnm={'fwyf1','mlrf1','lonf1','smkf1','sanf1','dryf1','looe1'};
  stnm=cstnm{:};
  if ( strcmpi(stnm,'looe1') )
    stn=optimize_station_heat_budget(stnm,'erai','avhrr_weekly','ndbc','tpxo_tide','erai',{'sfld','mc_seatemp'});
    stn=[]; clear stn;
    stn=optimize_station_heat_budget(stnm,'erai','avhrr_weekly','ndbc','tpxo_tide','erai',{'sfld','ad_seatemp'});
  else
    stn=optimize_station_heat_budget(stnm,'erai','avhrr_weekly','ndbc','tpxo_tide','erai');
  end;
  stn=[]; clear stn; pack;
end;

clear stnm cstnm

more on;
diary off;
