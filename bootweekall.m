function stats = bootweekall(stn)
%function stats = bootweekall(stn)
%
% Bootstrap estimate the weekly mean and stdev of sea and air temperature,
% spec. humid., wind speed, and three heat flux estimates: TOGA-COARE / NCEP,
% plus advection by Global HYCOM and Stokes drive, plus horizontal convection.
% Returns 3x52 mean and stdev. for each var., with 95% confidence interval.
%
% Last Saved Time-stamp: <Tue 2010-07-13 09:14:19  Lew.Gramer>

  datapath = get_thesis_path('../data');

  stnm = lower(stn.station_name);
  %DEBUG:
  disp(stnm);

  matfname = fullfile(datapath,sprintf('%s_bootweek.mat',stnm));
  if ( exist(matfname,'file') )

    load(matfname,'stats');

  else

    flds = { ...
        'ndbc_sea_t', ...
        'ndbc_air_t', ...
        'ndbc_relhumid', ...
        'ndbc_spechumid', ...
        'ndbc_wind1_speed', ...
        'ndbc_ncep_30_heat_flux_term', ...
        'ndbc_ncep_30_dt', ...
        'ndbc_ncep_30a_heat_flux_term', ...
        'ndbc_ncep_30a_dt', ...
        'ndbc_ncep_30a_absorbed_heat_flux', ...
        'ndbc_ncep_30a_absorbed_heat_flux_term', ...
        'ndbc_ncep_30a_total_heat_flux', ...
        'ndbc_ncep_30a_total_heat_flux_term', ...
        'netqf', ...
           };

    for fldix = 1:length(flds)
      fld = flds{fldix};
      %DEBUG:
      disp(fld);
      if ( ~isfield(stn,fld) || length(find(isfinite(stn.(fld).data))) == 0 )
        warning('No valid field "%s" in station struct.',fld);
      else
        tic,
        [stats.(fld).m,stats.(fld).s] = bootweek(stn,fld,300);
        toc,
      end;
    end;

    disp(['Saving to ' matfname]);
    save(matfname,'stats');

  end;

  plotbootweek(stats);

return;
