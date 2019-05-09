function dump_station_heat_budget_errors(stn,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,substitute_field_names)
%function dump_station_heat_budget_errors(stn,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,substitute_field_names)
%
% Output a report with DISP of the median and "variance" (max-min/6) of error
% in terms of heat budget as calculated by, e.g., STATION_HEAT_BUDGET_ERRORS.
%
% Last Saved Time-stamp: <Tue 2012-07-31 16:10:35  lew.gramer>

  if ( ~isfield(stn,'station_name') )
    error('Station struct STN had no .station_name field!');
  end;

  %%%
  %% Call SCRIPT to set:
  %% Set variable-name prefixes ("PFX") for various input/output datasets;
  %% AND, set all station struct fieldnames used to produce heat budget 
  %% NOTE WELL: Calls FIX_VARNAMELENGTHS to fix variable name string lengths
  %%             to meet MATLAB limitations. MAJOR SIDE EFFECTS...
  station_heat_budget_field_names;


  sQsh = stn.([qshfld '_err']).data(:);
  sQlh = stn.([qlhfld '_err']).data(:);
  sQlw = stn.([lrfld '_err']).data(:);
  sQsw = stn.([srfld '_err']).data(:);
  sgammaQsw = stn.([asrfld '_err']).data(:);
  sQ0 = stn.([sq0fld '_err']).data(:);
  sQt = stn.([sqtfld,'_err']).data(:);
  sgQ0 = stn.([q0fld '_err']).data(:);
  sgQ0t = stn.([q0fld '_err']).data(:);
  sQb = stn.([qbofld '_err']).data(:);
  sbQ0 = stn.([bq0fld '_err']).data(:);
  sbQ0t = stn.([bq0tfld,'_err']).data(:);

  sfqudT = stn.([fqudTfld,'_err']).data(:);
  skd2T = stn.([kd2Tfld,'_err']).data(:);
  sqtAdv = stn.([qtAdvfld,'_err']).data(:);
  sbdT = stn.([bdTfld,'_err']).data(:);
  shc = stn.([hcdTdt '_err']).data(:);

  sfqudTf = stn.([fqudTffld,'_err']).data(:);
  skd2Tf = stn.([kd2Tffld,'_err']).data(:);
  sqtAdvf = stn.([qtAdvffld,'_err']).data(:);
  sbdTf = stn.([bdTffld,'_err']).data(:);
  shcf = stn.([hcdTdtf,'_err']).data(:);


  disp('REPRESENTATION ERROR report');
  disp({'Qsh',nanmedian(sQsh),(nanmax(sQsh)-nanmin(sQsh))/6});
  disp({'Qlh',nanmedian(sQlh),(nanmax(sQlh)-nanmin(sQlh))/6});
  disp({'Qlw',nanmedian(sQlw),(nanmax(sQlw)-nanmin(sQlw))/6});
  disp({'Qsw',nanmedian(sQsw),(nanmax(sQsw)-nanmin(sQsw))/6});
  disp({'Q0',nanmedian(sQ0),(nanmax(sQ0)-nanmin(sQ0))/6});
  disp({'gQsw',nanmedian(sgammaQsw),(nanmax(sgammaQsw)-nanmin(sgammaQsw))/6});
  disp({'gQ0',nanmedian(sgQ0),(nanmax(sgQ0)-nanmin(sgQ0))/6});
  disp({'Qb',nanmedian(sQb),(nanmax(sQb)-nanmin(sQb))/6});
  disp({'gQ0+Qb',nanmedian(sbQ0),(nanmax(sbQ0)-nanmin(sbQ0))/6});
  disp({'fqudTf',nanmedian(sfqudTf),(nanmax(sfqudTf)-nanmin(sfqudTf))/6});
  disp({'kd2Tf',nanmedian(skd2Tf),(nanmax(skd2Tf)-nanmin(skd2Tf))/6});
  disp({'gQ0+Qb+fqudTf',nanmedian(sqtAdvf),(nanmax(sqtAdvf)-nanmin(sqtAdvf))/6});
  disp({'gQ0+Qb+fqudTf+kd2Tf',nanmedian(sbdTf),(nanmax(sbdTf)-nanmin(sbdTf))/6});
  disp({'dtTf',nanmedian(shcf),(nanmax(shcf)-nanmin(shcf))/6});
  disp('...');
  disp({'Q0/pCph',nanmedian(sQt),(nanmax(sQt)-nanmin(sQt))/6});
  disp({'Q0t=(gQ0+Qb)/pCph',nanmedian(sbQ0t),(nanmax(sbQ0t)-nanmin(sbQ0t))/6});
  disp({'fqudT',nanmedian(sfqudT),(nanmax(sfqudT)-nanmin(sfqudT))/6});
  disp({'kd2T',nanmedian(skd2T),(nanmax(skd2T)-nanmin(skd2T))/6});
  disp({'Q0t+fqudT',nanmedian(sqtAdv),(nanmax(sqtAdv)-nanmin(sqtAdv))/6});
  disp({'Q0t+fqudT+kd2T',nanmedian(sbdT),(nanmax(sbdT)-nanmin(sbdT))/6});
  disp({'dtT',nanmedian(shc),(nanmax(shc)-nanmin(shc))/6});

return;
