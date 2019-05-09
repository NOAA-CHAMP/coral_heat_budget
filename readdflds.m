function [station,cs] = readdflds(stnm,flds)
%function [station,cs] = readdflds(stnm,flds)

  datapath = get_thesis_path('../data');
  figspath = get_thesis_path('../figs');

  if ( ~exist('flds','var') || isempty(flds) )
    flds = { ...
        'nocs_net_heat_flux','nocs_heat_flux_term', ...
        'sat_net_heat_flux','sat_heat_flux_term','sat_wind_stress', ...
           };
  end;
  if ( ~iscellstr(flds) )
    flds = {flds};
  end;

  thfname = fullfile(datapath,sprintf('%s_thesis.mat',stnm)),
  if ( ~exist(thfname,'file') )
    error('File not found "%s"!',thfname);
  end;

  s = load(thfname);
  for fldix = 1:length(flds)
    fld = flds{fldix};
    if (isfield(s.(stnm),fld))
      x.(fld) = s.(stnm).(fld);
    else
      warning('Field "%s" not found in "%s"!',fld,thfname);
    end;
  end;
  s = []; clear s;


  dtfname = fullfile(datapath,sprintf('%s_dt.mat',stnm));
  if ( ~exist(dtfname,'file') )
    error('File not found "%s"!',dtfname);
  end;
  load(dtfname,'station');
  for fldix = 1:length(flds)
    fld = flds{fldix};
    if (isfield(x,fld))
      station.(fld) = x.(fld);
    end;
  end;
  x = []; clear x;

  switch ( stnm ),
   case 'lonf1',
    station = tryhc(station,true,[],0.15,0.35);
   case 'sanf1',  % STILL NOT SURE SETTINGS FOR THIS ONE??
    station = tryhc(station,true,[],1.20,0.35);
   otherwise,
    station = tryhc(station,false,[],0.30,0.30);
  end;

  msfname = fullfile(datapath,sprintf('%s_ms.mat',stnm));
  disp(['Saving ' msfname]);
  save(msfname,'station');

  cs = cmpclims(station);
  ylim('default');
  figfname = fullfile(figspath,sprintf('%s_cmpclims_netqf.png',stnm));
  print('-dpng',figfname);

return;
