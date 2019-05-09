function [grd,lons,lats] = grid_xyz(fname)

  grd = [];
  lons = [];
  lats = [];

  % Specs for data downloaded from 3-sec resolution NGDC Coastal Relief Model
  dl = 8.3300e-004;
  minlon = -85.0; maxlon = -78.0;
  minlat =  24.0; maxlat =  26.7;

  [grd, refvec] = spzerom([minlat maxlat], [minlon maxlon], (1/dl));

  % Load data from file, but only 10,000 rows/points at a time
  fid = fopen(fname,'r');
  while ( isempty(ferror(fid)) && ~feof(fid) )
    d = fscanf(fid,'%g,%g,%g\n',[3 10000])';
    badix = find(d(:,1) == 0 | d(:,2) == 0);
    d(badix,:) = [];
    if ( ~isempty(d) )
      lons = [lons d(:,1)'];
      lats = [lats d(:,2)'];
      grd = imbedm(d(:,2), d(:,1), d(:,3), grd, refvec);
    end;
    d = []; clear d;
  end;
  fclose(fid);

  lons = unique(lons);
  lats = unique(lats);

  % grd = cast(grd, 'single');

return;
