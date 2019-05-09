function lpath = get_thesis_path(ldir)
%function lpath = get_thesis_path(ldir)
%
% Store/retrieve data in this M-file's local directory

  %[pathroot, ig, ig, ig] = fileparts(mfilename('fullpath'));
  [pathroot, ig, ig] = fileparts(mfilename('fullpath'));

  if ( length(ldir) > 1 && strcmp(ldir(1:2), '..') )
    % This implements the "../" in a portable way...
    %[pathroot, ig, ig, ig] = fileparts(pathroot);
    [pathroot, ig, ig] = fileparts(pathroot);
    if ( length(ldir) == 2 )
      ldir = '';
    else
      ldir = ldir(3:end);
    end;
  end;

  lpath = fullfile(pathroot, ldir, '');

return;
