function testterms(cs,yr,begdy,enddy)
%function testterms(cs,yr,begdy,enddy)

  %linspec = {'b.','g.','r.','c.','k.','y.','m.','bd','gd','rd','cd','kd','yd','md'};
  %linspec = {'k.-','k.:','k.-.','k.--','ko-','ko:','ko-.','ko--','ks-','ks:','ks-.','ks--','k^-','k^:','k^-.','k^--'};
  %linspec = {'k.-','k.:','k.--','ko-','ko:','ko--','kp-','kp:','kp--','k^-','k^:','k^--','k*-','k*:','k*--'};
  linspec = {'k.-','k.:','ko-','ko:','kp-','kp:','k^-','k^:','k*-','k*:','kd-','kd:','ks-','ks:'};

  if ( nargin < 2 )
    yr=2005;
  end;
  if ( nargin < 3 )
    begdy = 1;
  end;
  if ( nargin < 4 )
    enddy = begdy+3;
  end;

  begdt = datenum(yr,1,1) + begdy - 1;
  enddt = datenum(yr,1,1) + enddy - 1;
  dts = begdt:(1/24):enddt-(29/(60*24));

  ys = [];
  for ix = 1:length(cs)
    csyr(ix) = find(cs(ix).yrs == yr);
    csbd(ix) = begdy*24;
    csed(ix) = (enddy*24) - 1;

    rawdat = cs(ix).datmtx(csyr(ix),csbd(ix):csed(ix));
    gdix = find(isfinite(rawdat));
    if ( ~isempty(gdix) )
      if ( gdix(1) ~= 1 )
        rawdat(1) = rawdat(gdix(1));
        gdix = [ 1 ; gdix(:) ];
      end;
      dat = interp1(dts(gdix),rawdat(gdix),dts,'linear',0);
    else
      dat = rawdat;
    end;

    if ( ix == 1 )
      dat = dat - dat(1);
    else
      dat = cumsum(dat);
    end;

    ys = [ ys ; dat ];
    legs(ix) = {cs(ix).name};
  end;

  legs = strrep(legs,'_','\_');

  figure;
  hold on;
  for cix = 1:length(cs)
    plot(dts,ys(cix,:),linspec{mod(cix-1,length(linspec))+1});
  end;
  maxigraph;
  set(gca,'position',[0.05 0.05 0.90 0.90]);
  datetick('x',0);
  set_datetick_cursor;
  legend(legs,'Location','Best');
  titlename( [datestr(begdt,1) ' - ' datestr(enddt,1) ]);
  grid on;

return;
