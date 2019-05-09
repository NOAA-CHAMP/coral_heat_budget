1;

diary(fullfile(get_thesis_path('../src'),'all_dump.log'));
timenow,
more off;

doPrint = false;

for cst={'fwyf1','mlrf1','lonf1','smkf1','sanf1','dryf1','looe1'};
  st=cst{:};
  disp(['%%%%%%%%%%%%%%%%%%%%']);
  disp(['%% Starting ',upper(st)]);
  disp(['%%%%%%%%%%%%%%%%%%%%']);
  if ( strcmpi(st,'looe1') )
    stn = optimize_station_heat_budget(st,'erai','avhrr_weekly','ndbc','tpxo_tide','erai',{'sfld','mc_seatemp'});
    if (doPrint)
      print('-dtiff',fullfile(get_thesis_path('../figs'),[lower(stn.station_name),'-mc-best.tif']));
    end;
    close all; stn=[]; clear stn;
    stn = optimize_station_heat_budget(st,'erai','avhrr_weekly','ndbc','tpxo_tide','erai',{'sfld','ad_seatemp'});
    if (doPrint)
      print('-dtiff',fullfile(get_thesis_path('../figs'),[lower(stn.station_name),'-ad-best.tif']));
    end;
    close all; stn=[]; clear stn;
  else
    stn = optimize_station_heat_budget(st,'erai','avhrr_weekly','ndbc','tpxo_tide','erai',{},struct());
    if (doPrint)
      print('-dtiff',fullfile(get_thesis_path('../figs'),[lower(stn.station_name),'-best.tif']));
    end;
    close all; stn=[]; clear stn;
  end;
end;

more on;
timenow,
diary off
