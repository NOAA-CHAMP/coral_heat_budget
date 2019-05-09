1;

more off;
%warning('off','anwera:NoURL');
warning('off','anwera:BadCoords');
warning('off','anwera:NoLocalAver');
warning('off','anwera:NoLocalAcc');
warning('off','anwera:NoData');
warning('off','anwera:BadData');

for yr=2005:2011
  switch (yr),
   case {2004,2008,2012,2016},
    jds = 1:366;
   case 2011,
    % Data service suspended
    jds = 1:33;
   otherwise,
    jds = 1:365;
  end;
  disp({yr,jds});
  for jd=jds(:)'
    for hr=0:23
      ds=sprintf('%04d%03d%02d00',yr,jd,hr);
      anwera(ds,0);
    end;
    pause(0.5);
  end;
end;

warning('on','anwera:NoURL');
warning('on','anwera:BadCoords');
warning('on','anwera:NoLocalAver');
warning('on','anwera:NoLocalAcc');
warning('on','anwera:NoData');
warning('on','anwera:BadData');
more on;
