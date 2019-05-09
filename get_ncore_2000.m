function stns = get_ncore_2000
%function stns = get_ncore_2000
%
% Load data for NCORE studies 2000-2002, provided courtesy of Dr. Tom N. Lee,
% Dr. Su Sponaugle, and Liz Williams of RSMAS. See Sponaugle et al 2003 and
% especially Sponaugle et al 2005.
%
% NOTE: Currently ignores "wandering" of C mooring deployment to deployment!
%
% Last Saved Time-stamp: <Thu 2013-03-28 16:50:49 Eastern Daylight Time Lew.Gramer>

 stns = [];

 matfname = fullfile(get_thesis_path('../data'),'ncore_2000.mat');
 if ( exist(matfname,'file') )
  disp(['Loading ',matfname]);
  load(matfname,'stns');

 else

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Load A/B/C moored current-meter/thermistor sites

  % Mooring "C" files actually from three successive locations: Stations
  % 'CSH' and 'CDP' combine all three locations into one set of time series.

  stnms = {'a','b','cdp','csh', 'cdp1','csh1','cdp2','csh2','cdp3','csh3'};

  % Cell array of cell-arrays of filenames: Order of the sub-arrays must
  % match the order of the cell array STNMS of station names above!
  st_fnames = { ...
      { ...
          'moorajun2000_3hlp.dat', ...
          'mooranov2000_3hlp.dat', ...
          'moora_part1_jun2001_3hlp.dat', ...
      }, ...
      { ...
          'moorbjun2000_3hlp.dat', ...
          'moorbjun2001_3hlp.dat', ...
      }, ...
      { ...
          'moorcbotjun2000_3hlp.dat', ...
          'moorcbotjun2001_3hlp.dat', ...
          'moorcbotoct2001_3hlp.dat', ...
          'moorcbotapr2002_3hlp.dat', ...
      }, ...
      { ...
          'moorctopjun2000_3hlp.dat', ...
          'moorctopjun2001_3hlp.dat', ...
          'moorctopoct2001_3hlp.dat', ...
          'moorctopapr2002_3hlp.dat', ...
      }, ...
      ...
      { ...
          'moorcbotjun2000_3hlp.dat', ...
      }, ...
      { ...
          'moorctopjun2000_3hlp.dat', ...
      }, ...
      { ...
          'moorcbotjun2001_3hlp.dat', ...
      }, ...
      { ...
          'moorctopjun2001_3hlp.dat', ...
      }, ...
      { ...
          'moorcbotoct2001_3hlp.dat', ...
          'moorcbotapr2002_3hlp.dat', ...
      }, ...
      { ...
          'moorctopoct2001_3hlp.dat', ...
          'moorctopapr2002_3hlp.dat', ...
      }, ...
  };


  % grep '\(Longitude\|Latitude\|Depth\)' moora*_3hlp.dat
  % moorajun2000_3hlp.dat:Latitude (+=N):                        25.1090
  % moorajun2000_3hlp.dat:Longitude (+=E):                      -80.3803
  % moorajun2000_3hlp.dat:Water Depth (m):                        4.1000
  % moorajun2000_3hlp.dat:Instrument Depth (m):                   0.0000
  % mooranov2000_3hlp.dat:Latitude (+=N):                        25.1090
  % mooranov2000_3hlp.dat:Longitude (+=E):                      -80.3803
  % mooranov2000_3hlp.dat:Water Depth (m):                        4.1000
  % mooranov2000_3hlp.dat:Instrument Depth (m):                   0.0000
  % moora_part1_jun2001_3hlp.dat:Latitude (+=N):                        25.1090
  % moora_part1_jun2001_3hlp.dat:Longitude (+=E):                      -80.3803
  % moora_part1_jun2001_3hlp.dat:Water Depth (m):                        4.1000
  % moora_part1_jun2001_3hlp.dat:Instrument Depth (m):                   0.0000

  % All files same location
  stns.a.lat = 25.1090;
  stns.a.lon = -80.3803;
  stns.a.depth = 4.1;
  stns.a.i_depth = 4.1;


  % grep '\(Longitude\|Latitude\|Depth\)' moorb*_3hlp.dat
  % moorbjun2000_3hlp.dat:Latitude (+=N):                        25.0932
  % moorbjun2000_3hlp.dat:Longitude (+=E):                      -80.3550
  % moorbjun2000_3hlp.dat:Water Depth (m):                        7.0000
  % moorbjun2000_3hlp.dat:Instrument Depth (m):                   0.0000
  %%%% moorbnov2000_3hlp.dat CONTAINED ALL -999.0 VALUES!
  %% moorbnov2000_3hlp.dat:Latitude (+=N):                        25.1090
  %% moorbnov2000_3hlp.dat:Longitude (+=E):                      -80.3803
  %% moorbnov2000_3hlp.dat:Water Depth (m):                        4.1000
  %% moorbnov2000_3hlp.dat:Instrument Depth (m):                   0.0000
  % moorbjun2001_3hlp.dat:Latitude (+=N):                        25.0932
  % moorbjun2001_3hlp.dat:Longitude (+=E):                      -80.3550
  % moorbjun2001_3hlp.dat:Water Depth (m):                        7.0000
  % moorbjun2001_3hlp.dat:Instrument Depth (m):                   0.0000

  % All files same location, except one with no valid data
  stns.b.lat = 25.0932;
  stns.b.lon = -80.3550;
  stns.b.depth = 7.0;
  stns.b.i_depth = 7.0;


  % grep '\(Longitude\|Latitude\|Depth\)' moorctop*_3hlp.dat
  % moorctopjun2000_3hlp.dat:Latitude (+=N):                        25.0740
  % moorctopjun2000_3hlp.dat:Longitude (+=E):                      -80.3178
  % moorctopjun2000_3hlp.dat:Water Depth (m):                       26.0000
  % moorctopjun2000_3hlp.dat:Instrument Depth (m):                   4.0000
  % moorctopjun2001_3hlp.dat:Latitude (+=N):                        25.0733
  % moorctopjun2001_3hlp.dat:Longitude (+=E):                      -80.3183
  % moorctopjun2001_3hlp.dat:Water Depth (m):                       21.8000
  % moorctopjun2001_3hlp.dat:Instrument Depth (m):                   0.0000
  % moorctopoct2001_3hlp.dat:Latitude (+=N):                        25.0673
  % moorctopoct2001_3hlp.dat:Longitude (+=E):                      -80.3183
  % moorctopoct2001_3hlp.dat:Water Depth (m):                       21.8000
  % moorctopoct2001_3hlp.dat:Instrument Depth (m):                   4.0000
  % moorctopapr2002_3hlp.dat:Latitude (+=N):                        25.0673
  % moorctopapr2002_3hlp.dat:Longitude (+=E):                      -80.3183
  % moorctopapr2002_3hlp.dat:Water Depth (m):                       21.8000
  % moorctopapr2002_3hlp.dat:Instrument Depth (m):                   4.0000

  % Files actually from three different locations: Stations 'CSH' and 'CDP'
  % combine all three locations into one set of time series.
  % 20000603 210000 to 20021002 150000 (with gaps)
  stns.csh.lat = 25.0673;
  stns.csh.lon = -80.3183;
  stns.csh.depth = 21.8;
  stns.csh.i_depth = 4.0;

  % 20000603 210000 to 20001109 120000
  stns.csh1.lat =  25.0740;
  stns.csh1.lon = -80.3178;
  stns.csh1.depth = 26.0;
  stns.csh1.i_depth = 4.0;
  % 20010620 000000 to 20011023 120000
  stns.csh2.lat =  25.0733;
  stns.csh2.lon = -80.3183;
  stns.csh2.depth = 21.8;
  stns.csh2.i_depth = 0.0;
  % 20011024 000000 to 20020404 150000 and 20020404 210000 to 20021002 150000
  stns.csh3.lat =  25.0673;
  stns.csh3.lon = -80.3183;
  stns.csh3.depth = 21.8;
  stns.csh3.i_depth = 4.0;


  % grep '\(Longitude\|Latitude\|Depth\)' moorcbot*_3hlp.dat
  % moorcbotjun2000_3hlp.dat:Latitude (+=N):                        25.0740
  % moorcbotjun2000_3hlp.dat:Longitude (+=E):                      -80.3178
  % moorcbotjun2000_3hlp.dat:Water Depth (m):                       26.0000
  % moorcbotjun2000_3hlp.dat:Instrument Depth (m):                  21.0000
  % moorcbotjun2001_3hlp.dat:Latitude (+=N):                        25.0733
  % moorcbotjun2001_3hlp.dat:Longitude (+=E):                      -80.3183
  % moorcbotjun2001_3hlp.dat:Water Depth (m):                       21.8000
  % moorcbotjun2001_3hlp.dat:Instrument Depth (m):                   0.0000
  % moorcbotoct2001_3hlp.dat:Latitude (+=N):                        25.0673
  % moorcbotoct2001_3hlp.dat:Longitude (+=E):                      -80.3183
  % moorcbotoct2001_3hlp.dat:Water Depth (m):                       21.8000
  % moorcbotoct2001_3hlp.dat:Instrument Depth (m):                  21.0000
  % moorcbotapr2002_3hlp.dat:Latitude (+=N):                        25.0673
  % moorcbotapr2002_3hlp.dat:Longitude (+=E):                      -80.3183
  % moorcbotapr2002_3hlp.dat:Water Depth (m):                       21.8000
  % moorcbotapr2002_3hlp.dat:Instrument Depth (m):                  21.0000

  % Files actually from three different locations: Stations 'CSH' and 'CDP'
  % combine all three locations into one set of time series.
  % 20000603 210000 to 20021002 150000 (with gaps)
  stns.cdp.lat = 25.0673;
  stns.cdp.lon = -80.3183;
  stns.cdp.depth = 21.8;
  stns.cdp.i_depth = 21.0;
  % 20000603 210000 to 20001109 090000
  stns.cdp1.lat =  25.0740;
  stns.cdp1.lon = -80.3178;
  stns.cdp1.depth = 26.0;
  stns.cdp1.i_depth = 21.0;
  % 20010620 000000 to 20011023 120000
  stns.cdp2.lat =  25.0733;
  stns.cdp2.lon = -80.3183;
  stns.cdp2.depth = 21.8;
  stns.cdp2.i_depth = 0;
  % 20011023 180000 to 20020404 150000 and 20020404 180000 to 20021010 210000
  stns.cdp3.lat =  25.0673;
  stns.cdp3.lon = -80.3183;
  stns.cdp3.depth = 21.8;
  stns.cdp3.i_depth = 21.0;


  % plot(-80.3178,25.0740,'wd'); text(-80.3178,25.0740,' <- jun2000','Color','w');
  % plot(-80.3183,25.0733,'wd'); text(-80.3183,25.0733,' <- jun2001','Color','w');
  % plot(-80.3183,25.0673,'wd'); text(-80.3183,25.0673,' <- oct2001*','Color','w');
  % plot(-80.3183,25.0673,'wd'); text(-80.3183,25.0673,' <- apr2002*','Color','w');


  for stix = 1:numel(stnms)

    stnm = stnms{stix};
    fnames = st_fnames{stix};

    t = []; u = []; v = []; T = [];

    for fnix = 1:numel(fnames)
      fname = fnames{fnix};
      fid = fopen(fullfile(get_thesis_path('../data'),'NCORE',fname),'r');
      if ( fid < 1 )
        warning('No file %s',fname);
      else
        for ln=1:22
          fgetl(fid);
        end;
        A = fscanf(fid,' %f %f %f %f %f\n');
        fclose(fid);

        n = numel(A)/5;
        if ( n ~= floor(n) )
          warning('Invalid format for file %s',fname);
        else
          dat = reshape(A',[5,n])';
          t(end+1:end+n) = datenum(dat(:,1),1,0) + dat(:,2);
          u(end+1:end+n) = dat(:,3);
          v(end+1:end+n) = dat(:,4);
          T(end+1:end+n) = dat(:,5);
        end; %if n else
      end; %if fid else
    end; %for fnix

    stns.(stnm).cm_u_3hlp.date = t(:);
    stns.(stnm).cm_u_3hlp.data = u(:);
    stns.(stnm).cm_v_3hlp.date = t(:);
    stns.(stnm).cm_v_3hlp.data = v(:);
    stns.(stnm).cm_seatemp_3hlp.date = t(:);
    stns.(stnm).cm_seatemp_3hlp.data = T(:);

    % 42 degT chosen as best guess for alongshore direction at all sites
    stns.(stnm) = station_reorient_vectors(stns.(stnm),42,'cm_u_3hlp','cm_v_3hlp');

  end; %for stix



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Load TL* bottom-mounted thermistor sites

  % From $THESISPATH/data/NCORE/tl_deploy1/ncore1tab.doc
  % TL1	25 04.14	80 19.17	9.75	TKSA	9.75	050	5	06/03/2000 1923	11/09/2000 1420	no data
  % TL2	25 04.43	80 19.47	7.0	TKSA	7.0	056	5	06/03/2000 1901	11/09/2000 1406	good
  % TL3	25 04.77	80 20.06	4.0	TKSA	4.0	046	5	06/03/2000 1847	11/08/2000 2056	no data
  % TL4	25 05.27	80 20.83	5.6	TKSA	5.6	049	5	06/03/2000 1830	11/08/2000 2030	good
  % From $THESISPATH/data/NCORE/tl_deploy2/ncore2tab.doc
  % TL1	25 04.14	80 19.17	9.75	TKSA	9.75	TL112	5	11/09/2000 1417	06/19/2001 1514	good
  % TL2	25 04.43	80 19.47	7.0	TKSA	7.0	TL23	5	11/09/2000 1402	06/19/2001 1527	good
  % TL3	25 04.77	80 20.06	4.0	TKSA	4.0	TL92	5	11/08/2000 2049	06/19/2001 1539	good
  % TL4	25 05.27	80 20.83	5.6	TKSA	5.6	TL17	5	11/08/2000 2024	06/19/2001 1551	good
  % From $THESISPATH/data/NCORE/tl_deploy3/ncore3tab.doc
  % TL1	25 04.14	80 19.17	9.75	TKSA	9.75	TL20	5	06/19/2001 1514	10/23/2001 1801	good
  % TL2	25 04.43	80 19.47	7.0	TKSA	7.0	Tl46	5	06/19/2001 1527	10/23/2001 1813	good
  % TL3	25 04.77	80 20.06	4.0	TKSA	4.0	TL49	5	06/19/2001 1539	10/23/2001 1821	good
  % TL4	25 05.27	80 20.83	5.6	TKSA	5.6	TL56	5	06/19/2001 1551	10/23/2001 1829	good
  % From $THESISPATH/data/NCORE/tl_deploy4/ncore4tab.doc
  % TL1	25 04.14	80 19.17	9.75	TKSA	9.75	TL103	5	10/23/2001 1801	07/10/2002 1706	
  % TL2	25 04.43	80 19.47	7.0	TKSA	7.0	TL92	5	10/23/2001 1813	07/10/2002 1725	
  % TL3	25 04.77	80 20.06	4.0	TKSA	4.0	TL23	5	10/23/2001 1821	07/10/2002 1744	
  % TL4	25 05.27	80 20.83	5.6	TKSA	5.6	TL26	5	10/23/2001 1829	07/10/2002 1800	
  % From $THESISPATH/data/NCORE/tl_deploy5/ncore5tab.doc
  % TL1	25 04.14	80 19.17	9.75	TKSA	9.75	TL49	5	07/10/2002 1654	11/21/2002 1618	
  % TL2	25 04.43	80 19.47	7.0	TKSA	7.0	TL83	5	07/10/2002 1712	11/21/2002 1628	
  % TL3	25 04.77	80 20.06	4.0	TKSA	4.0	TL82	5	07/10/2002 1736	11/21/2002 1641	
  % TL4	25 05.27	80 20.83	5.6	TKSA	5.6	TL57	5	07/10/2002 1755	11/20/2002 1923	

  stnms = {'tl1','tl2','tl3','tl4'};
  st_fnames = { ...
      { ...
          'TL1_2.txt', ...
          'TL1_3.txt', ...
          'TL1_5.txt', ...
      }, ...
      { ...
          'TL2_1.txt', ...
          'TL2_2.txt', ...
          'TL2_3.txt', ...
          'TL2_4.txt', ...
          'TL2_5.txt', ...
      }, ...
      { ...
          'TL3_2.txt', ...
          'TL3_3.txt', ...
          'TL3_4.txt', ...
          'TL3_5.txt', ...
      }, ...
      { ...
          'TL4_1.txt', ...
          'TL4_2.txt', ...
          'TL4_3.txt', ...
          'TL4_4.txt', ...
          'TL4_5.txt', ...
      }, ...
  };


  % grep -m 1 '  [0-9][0-9][.]' TL*_*.txt
  % TL1_2.txt:16:25:10       11/3/00
  % TL1_3.txt:12:30:16       6/15/01
  % TL1_5.txt:10:00:03       6/21/02 1st REC
  % TL2_1.txt:12:30:11       6/1/00
  % TL2_2.txt:14:45:10       11/3/00
  % TL2_3.txt:12:50:22       6/15/01
  % TL2_4.txt:12:40:03       10/8/01
  % TL2_5.txt:12:40:11       6/21/02 
  % TL3_2.txt:16:10:10       11/3/00
  % TL3_3.txt:13:30:17       6/15/01
  % TL3_4.txt:13:40:04       10/8/01
  % TL3_5.txt:12:15:13       6/21/02         
  % TL4_1.txt:13:30:08       6/1/00
  % TL4_2.txt:17:00:09       11/3/00
  % TL4_3.txt:13:40:12       6/15/01
  % TL4_4.txt:13:20:04      
  % Initial date for TL4_4.txt is a guess based on initial dates for TL[23]_4.txt
  % TL4_4.txt:13:20:04       10/8/01
  % TL4_5.txt:11:30:12       6/21/02 1st REC
  dtfmt = 'HH:MM:SS m/d/yy';
  st_inidts = { ...
      { ...
          datenum('16:25:10       11/3/00',dtfmt), ... %TL1_2.txt
          datenum('12:30:16       6/15/01',dtfmt), ... %TL1_3.txt
          datenum('10:00:03       6/21/02',dtfmt), ... %TL1_5.txt
      }, ...
      { ...
          datenum('12:30:11       6/1/00  ',dtfmt), ... %TL2_1.txt
          datenum('14:45:10       11/3/00 ',dtfmt), ... %TL2_2.txt
          datenum('12:50:22       6/15/01 ',dtfmt), ... %TL2_3.txt
          datenum('12:40:03       10/8/01 ',dtfmt), ... %TL2_4.txt
          datenum('12:40:11       6/21/02 ',dtfmt), ... %TL2_5.txt
      }, ...
      { ...
          datenum('16:10:10       11/3/00         ',dtfmt), ... %TL3_2.txt
          datenum('13:30:17       6/15/01         ',dtfmt), ... %TL3_3.txt
          datenum('13:40:04       10/8/01         ',dtfmt), ... %TL3_4.txt
          datenum('12:15:13       6/21/02         ',dtfmt), ... %TL3_5.txt
      }, ...
      { ...
          datenum('13:30:08       6/1/00 ',dtfmt), ... %TL4_1.txt
          datenum('17:00:09       11/3/00',dtfmt), ... %TL4_2.txt
          datenum('13:40:12       6/15/01',dtfmt), ... %TL4_3.txt
          datenum('13:20:04       10/8/01',dtfmt), ... %TL4_4.txt
          datenum('11:30:12       6/21/02',dtfmt), ... %TL4_5.txt
      }, ...
  };

  % Dates to subset for quality control
  dtfmt = 'mm/dd/yyyy HHMM';
  % TL1	25 04.14	80 19.17	9.75	TKSA	9.75	TL112	5	11/09/2000 1417	06/19/2001 1514	good
  % TL1	25 04.14	80 19.17	9.75	TKSA	9.75	TL20	5	06/19/2001 1514	10/23/2001 1801	good
  % TL1	25 04.14	80 19.17	9.75	TKSA	9.75	TL49	5	07/10/2002 1654	11/21/2002 1618	

  % TL2	25 04.43	80 19.47	7.0	TKSA	7.0	056	5	06/03/2000 1901	11/09/2000 1406	good
  % TL2	25 04.43	80 19.47	7.0	TKSA	7.0	TL23	5	11/09/2000 1402	06/19/2001 1527	good
  % TL2	25 04.43	80 19.47	7.0	TKSA	7.0	Tl46	5	06/19/2001 1527	10/23/2001 1813	good
  % TL2	25 04.43	80 19.47	7.0	TKSA	7.0	TL92	5	10/23/2001 1813	07/10/2002 1725	
  % TL2	25 04.43	80 19.47	7.0	TKSA	7.0	TL83	5	07/10/2002 1712	11/21/2002 1628	

  % TL3	25 04.77	80 20.06	4.0	TKSA	4.0	TL92	5	11/08/2000 2049	06/19/2001 1539	good
  % TL3	25 04.77	80 20.06	4.0	TKSA	4.0	TL49	5	06/19/2001 1539	10/23/2001 1821	good
  % TL3	25 04.77	80 20.06	4.0	TKSA	4.0	TL23	5	10/23/2001 1821	07/10/2002 1744	
  % TL3	25 04.77	80 20.06	4.0	TKSA	4.0	TL82	5	07/10/2002 1736	11/21/2002 1641	

  % TL4	25 05.27	80 20.83	5.6	TKSA	5.6	049	5	06/03/2000 1830	11/08/2000 2030	good
  % TL4	25 05.27	80 20.83	5.6	TKSA	5.6	TL17	5	11/08/2000 2024	06/19/2001 1551	good
  % TL4	25 05.27	80 20.83	5.6	TKSA	5.6	TL56	5	06/19/2001 1551	10/23/2001 1829	good
  % TL4	25 05.27	80 20.83	5.6	TKSA	5.6	TL26	5	10/23/2001 1829	07/10/2002 1800	
  % TL4	25 05.27	80 20.83	5.6	TKSA	5.6	TL57	5	07/10/2002 1755	11/20/2002 1923	
  st_begdts = { ...
      { ...
          datenum('11/09/2000 1417',dtfmt), ... %TL1_2.txt
          datenum('06/19/2001 1514',dtfmt), ... %TL1_3.txt
          datenum('07/10/2002 1654',dtfmt), ... %TL1_5.txt
      }, ...
      { ...
          datenum('06/03/2000 1901',dtfmt), ... %TL2_1.txt
          datenum('11/09/2000 1402',dtfmt), ... %TL2_2.txt
          datenum('06/19/2001 1527',dtfmt), ... %TL2_3.txt
          datenum('10/23/2001 1813',dtfmt), ... %TL2_4.txt
          datenum('07/10/2002 1712',dtfmt), ... %TL2_5.txt
      }, ...
      { ...
          datenum('11/08/2000 2049',dtfmt), ... %TL3_2.txt
          datenum('06/19/2001 1539',dtfmt), ... %TL3_3.txt
          datenum('10/23/2001 1821',dtfmt), ... %TL3_4.txt
          datenum('07/10/2002 1736',dtfmt), ... %TL3_5.txt
      }, ...
      { ...
          datenum('06/03/2000 1830',dtfmt), ... %TL4_1.txt
          datenum('11/08/2000 2024',dtfmt), ... %TL4_2.txt
          datenum('06/19/2001 1551',dtfmt), ... %TL4_3.txt
          datenum('10/23/2001 1829',dtfmt), ... %TL4_4.txt
          datenum('07/10/2002 1755',dtfmt), ... %TL4_5.txt
      }, ...
  };
  st_enddts = { ...
      { ...
          datenum('06/19/2001 1514',dtfmt), ... %TL1_2.txt
          datenum('10/23/2001 1801',dtfmt), ... %TL1_3.txt
          datenum('11/21/2002 1618',dtfmt), ... %TL1_5.txt
      }, ...
      { ...
          datenum('11/09/2000 1406',dtfmt), ... %TL2_1.txt
          datenum('06/19/2001 1527',dtfmt), ... %TL2_2.txt
          datenum('10/23/2001 1813',dtfmt), ... %TL2_3.txt
          datenum('07/10/2002 1725',dtfmt), ... %TL2_4.txt
          datenum('11/21/2002 1628',dtfmt), ... %TL2_5.txt
      }, ...
      { ...
          datenum('06/19/2001 1539',dtfmt), ... %TL3_2.txt
          datenum('10/23/2001 1821',dtfmt), ... %TL3_3.txt
          datenum('07/10/2002 1744',dtfmt), ... %TL3_4.txt
          datenum('11/21/2002 1641',dtfmt), ... %TL3_5.txt
      }, ...
      { ...
          datenum('11/08/2000 2030',dtfmt), ... %TL4_1.txt
          datenum('06/19/2001 1551',dtfmt), ... %TL4_2.txt
          datenum('10/23/2001 1829',dtfmt), ... %TL4_3.txt
          datenum('07/10/2002 1800',dtfmt), ... %TL4_4.txt
          datenum('11/20/2002 1923',dtfmt), ... %TL4_5.txt
      }, ...
  };

  stns.tl1.lat =   25.06900;
  stns.tl1.lon =  -80.31950;
  stns.tl1.depth = 9.75;
  stns.tl1.i_depth = 9.75;

  stns.tl2.lat =   25.07383;
  stns.tl2.lon =  -80.32450;
  stns.tl2.depth = 7.0;
  stns.tl2.i_depth = 7.0;

  stns.tl3.lat =   25.07950;
  stns.tl3.lon =  -80.33433;
  stns.tl3.depth = 4.0;
  stns.tl3.i_depth = 4.0;

  stns.tl4.lat =   25.08783;
  stns.tl4.lon =  -80.34717;
  stns.tl4.depth = 5.6;
  stns.tl4.i_depth = 5.6;

  for stix = 1:numel(stnms)
    stnm = stnms{stix};
    fnames = st_fnames{stix};
    inidts = st_inidts{stix};
    begdts = st_begdts{stix};
    enddts = st_enddts{stix};

    t = []; T = [];

    for fnix = 1:numel(fnames)
      fname = fnames{fnix};
      inidt = inidts{fnix};
      begdt = begdts{fnix};
      enddt = enddts{fnix};

      fid = fopen(fullfile(get_thesis_path('../data'),'NCORE',fname),'r');
      if ( fid < 1 )
        warning('No file %s',fname);
      else
        % With cross-file variations, processing this format is very ugly...
        A = textscan(fid,'%*[^\t] %f%*s');
        fclose(fid);

        dat = A{1};
        n = numel(dat);
        if ( n < 100 )
          warning('Invalid format for file %s',fname);
        else
          % Sample interval for all TL deployments was 5 mins.
          dts = inidt + ([0:n-1]*(5/60/24));

          % Quality control
          goodix = find(begdt<=dts & dts <=enddt);
          dat = dat(goodix(1+4:end-4));
          dts = dts(goodix(1+4:end-4));
          n = numel(dat);

          t(end+1:end+n) = dts;
          T(end+1:end+n) = dat;
        end; %if n else
        dts=[]; dat=[]; A=[]; clear dts dat A
      end; %if fid else
    end; %for fnix

    stns.(stnm).tl_seatemp.date = t(:);
    stns.(stnm).tl_seatemp.data = T(:);

    % 36 5-minute datapoints per 3 hours
    [stns.(stnm).tl_seatemp_3hlp.date,stns.(stnm).tl_seatemp_3hlp.data] = lanczos_ts(t(:),T(:),36);
  end; %for stix


  % Save to a MAT file for ease of future reference
  disp(['Saving ',matfname]);
  save(matfname,'stns');

 end; %if ( exist(matfname,'file') ) else

return;
