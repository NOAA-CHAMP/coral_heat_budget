function togo(tilwhat)
%function togo(tilwhat)
% Time left in PhD

  if ( ~exist('tilwhat','var') || isempty(tilwhat) )
    tilwhat = 'draft';
  end;
  switch (lower(tilwhat)),
   case 'kml',			tgt=datenum(2013, 9,01, 9, 0, 0);
   case 'crcp',			tgt=datenum(2013, 8,15,16, 0, 0);
   case 'done',			tgt=datenum(2013, 8,02,16, 0, 0);
   case 'final',		tgt=datenum(2013, 7,24,16, 0, 0);
   case 'johnt',		tgt=datenum(2013, 7,18,11, 0, 0);
   case 'defense',		tgt=datenum(2013, 7,16,11, 0, 0);
   case 'announcement',	tgt=datenum(2013, 7, 2,11, 0, 0);
   case 'draft',		tgt=datenum(2013, 6,24,17, 0, 0); disp('Do not forget Announcement form!');
   case 'trip',			tgt=datenum(2013, 6,13,20, 0, 0);
   otherwise,			tgt=datenum(2013, 6,24,17, 0, 0);
  end;
  x=tgt-now;
  disp(sprintf('%gd %gh',floor(x),floor((x-floor(x))*24)));
return;
