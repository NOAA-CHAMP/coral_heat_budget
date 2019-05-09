function stats = bootweekadd(stn)
%function stats = bootweekadd(stn)
%
% Bootstrap estimate the weekly mean and stdev of sea and air temperature,
% spec. humid., wind speed, and three heat flux estimates: TOGA-COARE / NCEP,
% plus advection by Global HYCOM and Stokes drive, plus horizontal convection.
% Returns 3x52 mean and stdev. for each var., with 95% confidence interval.
%
% Last Saved Time-stamp: <Tue 2010-07-13 09:35:58  Lew.Gramer>

  datapath = get_thesis_path('../data');

  stnm = lower(stn.station_name);
  %DEBUG:
  disp(stnm);

  matfname = fullfile(datapath,sprintf('%s_bootweek.mat',stnm));

  load(matfname,'stats');

      % 'ndbc_ncep_30a_absorbed_heat_flux', ...
      % 'ndbc_ncep_30a_absorbed_heat_flux_term', ...
      % 'ndbc_ncep_30a_total_heat_flux', ...
      % 'ndbc_ncep_30a_total_heat_flux_term', ...
  flds = { ...
      'netqf', ...
         };

  for fldix = 1:length(flds)
    fld = flds{fldix};
    if ( isfield(stn,fld) )
      %DEBUG:
      disp(fld);
      tic,
      % saveds = stats.(fld).s;
      [stats.(fld).m,stats.(fld).s] = bootweek(stn,fld,300);
      % stats.(fld).s([1 3],:) = saveds([1 3],:);
      toc,
    end;
  end;

  disp(['Resaving to ' matfname]);
  save(matfname,'stats');

  plotbootweek(stats);

return;
