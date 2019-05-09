function result = read_all_insolation_csvs(stanm)
%function result = read_all_insolation_csvs(stanm)
%
% Read all available CSV files with up/down short- (insolation) and longwave
% radiation, PAR, and other values estimated by NESDIS from satellite data,
% for gridpoints surrounding the reef monitoring site named 'STANM'. RESULT
% is usually a struct with fields .DATE, .PAR, .SU, .SD, .LU, .LD, etc.
%
% Last Saved Time-stamp: <Mon 2011-10-24 16:00:54  lew.gramer>

  result = [];

  datapath = get_thesis_path('../data');

  for yr = 1993:get_year(now)

    fname = fullfile(datapath, sprintf('insolation_%s_%d.csv', stanm, yr));

    if ( exist(fname,'file') )

      disp(yr);
      dat = read_insolation_csv(fname, stanm);
      flds = fieldnames(dat);
      for ix = 1:length(flds)
        fld = flds{ix};
        if ( ~isfield(result, fld) )
          result.(fld) = [];
        end;
        result.(fld)(end+1:end+numel(dat.(fld))) = dat.(fld);
      end;

    end;

  end;

  if ( isempty(result) )
    warning('Found NO insolation CSV files, e.g., no "%s"!', fname);
  end;

return;
