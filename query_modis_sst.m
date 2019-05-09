function [LON,LAT,VAL] = query_modis_sst(bbox)

error('UNDER CONSTRUCTION');

  if ( historical_data )
    querysst = sprintf('%s.ascii?sst%s', ...
                       sstfname, ixstr);
  else
    querysst = sprintf('%s.ascii?seadas_sst%s,l2_flags%s,cloud_mask%s', ...
                       sstfname, ixstr, ixstr, ixstr);
  end;
  % DEBUG:  disp('SST:'); toc;
  sst_s = query_dods(baseurl, querysst, nrows, ncols);
  if ( isempty(sst_s) )
    querysst = sprintf('%s.ascii?seadas_sst%s,cloud_mask%s', ...
                       sstfname, ixstr, ixstr);
    sst_s = query_dods(baseurl, querysst, nrows, ncols);
    if ( isempty(sst_s) )
      querysst = sprintf('%s.ascii?seadas_sst%s', ...
                         sstfname, ixstr);
      sst_s = query_dods(baseurl, querysst, nrows, ncols);
    end;
  end;
  % DEBUG:  disp('Done'); toc;


return;
