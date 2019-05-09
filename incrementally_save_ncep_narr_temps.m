1;

more off;

datapath = get_thesis_path('../data');

ALLVARS = { ...
    'Temperature_surface', ...
    'Temperature_height_above_ground', ...
          };

%for yr = [ 2010 2009 2008 2007 2006 2005 2004 2003 2002 2001 2000 1999 1998 1997 1996 1995 1994 1993 1992 1991 1990 1989 1988 1987 ]
for yr = [ ]

 switch ( yr )
  case 2010,
   mos = 1:3;
  otherwise,
   mos = 1:12;
 end;

 for mo = mos(:)'

  dat = []; clear dat;

  fname = fullfile(datapath, sprintf('ncep_narr_temps_%04d_%02d.mat', yr, mo));

  disp(['Creating ' fname '...']);

  dat = query_ncep_narr_subset([],datenum(yr,mo,1),(datenum(yr,(mo+1),1)-1),ALLVARS);
  disp(['Saving dat to ' fname]);
  save(fname, 'dat');

  % Subset our world list of stations to those inside our BBOX
  if ( ~exist('stns','var') || isempty(stns) )
    stns = get_all_station_metadata;
    XV = [dat.bbox(3) dat.bbox(4) dat.bbox(4) dat.bbox(3)];
    YV = [dat.bbox(2) dat.bbox(2) dat.bbox(1) dat.bbox(1)];
    goodix = find( inside(stns.lons, stns.lats, XV, YV) );
    clear XV YV;
    stns.codes  = stns.codes(goodix);
    stns.lons   = stns.lons(goodix);
    stns.lats   = stns.lats(goodix);
    stns.depths = stns.depths(goodix);
  end;


  for ix = 1:length(stns.lons)

    stnm = lower(stns.codes{ix});
    disp(stnm);
    [lonix,latix] = gridnbhd_km(dat.lon,dat.lat,stns.lons(ix),stns.lats(ix),0);
    dat.(stnm).lonix = lonix;
    dat.(stnm).latix = latix;

    dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'sea_t', 'Temperature_surface');
    dat.(stnm) = extract_ncep_field(dat, dat.(stnm), 'air_t', 'Temperature_height_above_ground');

  end;

  %%%% ??? DEBUG
  disp(['Saving again to ' fname]);
  save(fname, 'dat');

 end; %for mo

end; %for yr
