function stats = bootmonadd(stn)
%function stats = bootmonadd(stn)
%
% Bootstrap estimate the monthly mean and stdev of sea and air temperature,
% spec. humid., wind speed, and three heat flux estimates: TOGA-COARE / NCEP,
% plus advection by Global HYCOM and Stokes drive, plus horizontal convection.
% Returns 3x12 mean and stdev. for each var., with 95% confidence interval.
%
% Last Saved Time-stamp: <Tue 2010-07-13 08:12:12  Lew.Gramer>

  datapath = get_thesis_path('../data');

  stnm = lower(stn.station_name);
  %DEBUG:
  disp(stnm);

  matfname = fullfile(datapath,sprintf('%s_bootmon.mat',stnm));

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
      [stats.(fld).m,stats.(fld).s] = bootmon(stn,fld,300);
    end;
  end;

  disp(['Resaving to ' matfname]);
  save(matfname,'stats');

  plotbootmon(stats);

return;
