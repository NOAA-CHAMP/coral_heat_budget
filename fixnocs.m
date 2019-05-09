function fixnocs
%  error('This was a RUN-ONCE fixit function: you should not run this!');
error('This was a RUN-ONCE fixit function: you should not run this!');

  datapath = get_thesis_path('../data');

  for stnmc = { 'lkwf1','fwyf1','mlrf1','lonf1','tnrf1','smkf1','sanf1','dryf1' }
    stnm = stnmc{:};
    matfname = fullfile(datapath,sprintf('%s_nocs.mat',stnm));
    %DEBUG:
    disp(matfname);
    load(matfname,'stn');
    result = stn;
    stn = []; clear stn;
    delete(matfname);
    save(matfname,'result');
  end;

return;
