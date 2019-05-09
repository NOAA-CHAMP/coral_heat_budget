1;

set_more off;

diary get_erai_station.log

stnms = { 'bnpnn','bnpon','bnppa','bnpmi','bnpin','condp','consh','tavrk',...
          'lciy2','42003','dryf1','plsf1','sanf1','amsf1','looe1','smkf1',...
          'mose1','tnrf1','lonf1','mlrf1','cryf1','fwyf1','lkwf1','41140',...
          'dbjm1','lppr1','srvi2','cmrc3' };

stns = get_erai_station(stnms);

% for cstnm={'cryf1','dryf1','fwyf1','lkwf1','lonf1','looe1','mlrf1','mose1','plsf1','sanf1','smkf1','tnrf1'};
%   stnm=cstnm{:};
%   stn = get_station_from_station_name(stnm);
%   stn = get_ww3_station(stn);
%   stn=[]; clear stn;
% end;

diary off

set_more;
