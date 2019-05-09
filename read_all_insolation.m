function stn = read_all_insolation(stn)
%function stn = read_all_insolation(stn)
%
% Read all available CSV files with up/down short- (insolation) and longwave
% radiation, PAR, and other values estimated by NESDIS from satellite data,
% for gridpoints surrounding reef monitoring site named 'STN.station_name'.
% STN is returned with fields DATE, PAR, SU, SD, LU, LD, etc.
%
% Last Saved Time-stamp: <Mon 2010-03-22 12:12:55 Eastern Daylight Time Lew.Gramer>

  stanm = stn.station_name;

  result = read_all_insolation_csvs(stanm);

  gsippath = get_thesis_path('../data/class_gsip');

  % fknms_20090402_0045.v7.3.mat
  matfiles = dir(fullfile(gsippath,'fknms*.v7.3.mat'));
  for fidx = 1:length(matfiles)
    matfname = fullfile(gsippath, [matfiles(fidx).name]);
    %disp(sprintf('Loading MAT sat insolation file %s', matfname));
    load(matfname);
    if ( ~exist('lonix','var') )
      lonix = find(dat.lon
    end;
  end;

return;
