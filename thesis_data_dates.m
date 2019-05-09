1;

more off;

% Retrieve data from a path relative to this M-file's local directory
if ( ~exist('datapath', 'var') || isempty(datapath) )
  datapath = get_thesis_path('../data');
end;

stnms = {'LKWF1','FWYF1','MLRF1','LONF1','SMKF1','SANF1','DRYF1','42003'};

for ix = 1:length(stnms)
  stnm = stnms{ix};
  dts{ix}.station_name = stnm;

  station = []; clear station;
  % station = load_all_ndbc_data([], stnm);
  % save(fullfile(datapath, [stnm '-ndbc.mat']), 'station');
  disp(['Loading ' stnm '-ndbc.mat']);
  load(fullfile(datapath, [stnm '-ndbc.mat']), 'station');

  dts{ix}.wyrs = [];  dts{ix}.tyrs = [];  dts{ix}.dyrs = [];  dts{ix}.pyrs = [];

  if ( isfield(station, 'wind1_speed') );
    dvec = datevec(station.wind1_speed.date);
    dts{ix}.wyrs = unique(dvec(:,1));
  end;
  if ( isfield(station, 'air_t') );
    dvec = datevec(station.air_t.date);
    dts{ix}.ayrs = unique(dvec(:,1));
  end;
  if ( isfield(station, 'sea_t') );
    dvec = datevec(station.sea_t.date);
    dts{ix}.tyrs = unique(dvec(:,1));
  end;
  if ( isfield(station, 'dew_t') );
    dvec = datevec(station.dew_t.date);
    dts{ix}.dyrs = unique(dvec(:,1));
  end;
  if ( isfield(station, 'licor_surf_par') );
    dvec = datevec(station.licor_surf_par.date);
    dts{ix}.pyrs = unique(dvec(:,1));
  end;
end;

station = []; clear station;
clear datapath stnms stnm ix dvec;

fprintf( 1, '\n\nREPORT:\n' );
for ix = 1:length(dts)
  fprintf( 1, '\n %s: W:%d-%d A:%d-%d T:%d-%d D:%d-%d P:%d-%d\n', ...
           dts{ix}.station_name, ...
           min(dts{ix}.wyrs), max(dts{ix}.wyrs), ...
           min(dts{ix}.ayrs), max(dts{ix}.ayrs), ...
           min(dts{ix}.tyrs), max(dts{ix}.tyrs), ...
           min(dts{ix}.dyrs), max(dts{ix}.dyrs), ...
           min(dts{ix}.pyrs), max(dts{ix}.pyrs) );
end;

more on;
