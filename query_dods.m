function s = query_dods(baseurl, queryargs, nrows, ncols)
%function s = query_dods(baseurl, queryargs, nrows, ncols)
%
% Submit a query to a DoDS (text access to binary data subset) server. Second
% arg 'queryargs' should have no file separators ("/" or "\") in it, so that
% it is suitable for conversion into a local filename where the ASCII query
% result will be cached for future requests. Returns a struct 's' containing
% one field for each variable succesfully queried by the 'queryargs' string.
% Each field will contain a numeric array - of double for most variables, but
% of 'uint32' if the variable name happens to be 'l2_flags' (HACK ALERT).
%
% Last Saved Time-stamp: <Mon 2010-10-18 11:38:45 Eastern Daylight Time gramer>

  s = [];

  datapath = get_thesis_path('../data');
  modispath = fullfile(datapath,'modis');
  queryfname = fullfile(modispath, [regexprep(queryargs, '\W+', '_') '.txt']);

  % DEBUG:  fprintf(1, 'querystr=%s, ncols=%d\n', querystr, ncols);
  % DEBUG:  tic;
  if ( ~exist(queryfname, 'file') )
    querystr = sprintf('%s/%s', baseurl, queryargs);
    % DEBUG:
    disp(['Trying to download ' querystr]);
    [f, status] = urlwrite(querystr, queryfname);
    if ( isempty(f) || status ~= 1 )
      warning('URLWRITE("%s") failed!', querystr);
      return;
    end;
  end;
  % DEBUG:  disp('urlwrite'); toc;

  fmt = repmat(',%f', [1 ncols]); fmt = [ '%[^[][%d]' fmt ];
  flgfmt = repmat(',%d32', [1 ncols]); flgfmt = [ '%[^[][%d]' flgfmt ];

  % Parse each line of result file separately
  fid = fopen(queryfname, 'r');
  if ( fid < 0 )
    error('Cannot open cached query-result file %s', queryfname);
  end;
  res = textscan(fid, '%[^\n]\n');
  fclose(fid);
  % DEBUG:  disp('textscan'); toc;

  res = res{1};
  for nline = 1:length(res)
    ln = res{nline};
    % DEBUG:    if (mod(nline,100)==0); disp('start'); toc; end;
    if ( ~isempty(ln) )
      clear cstrs;
      cstrs = textscan(ln, fmt);
      % DEBUG:      if (mod(nline,100)==0); disp('textscan'); toc; end;
      if ( length(cstrs) >= 3 && ~isempty(cstrs{3}) )
        var = cstrs{1}{:};
        % row = cstrs{2} + 1;
        row = nrows - cstrs{2};

        % Handle bit flags specially
        % (Idiotic mazzafrazzin MATLAB)
        if ( strcmpi(var, 'l2_flags') )
          if ( ~isfield(s, var) )
            s.(var) = repmat(intmax('uint32'), [nrows ncols]);
          end;
          clear cstrs;
          cstrs = textscan(ln, flgfmt);
          s.(var)(row,1:length(cstrs)-2) = typecast([cstrs{3:end}],'uint32');

        else
          if ( ~isfield(s, var) )
            s.(var) = repmat(nan, [nrows ncols]);
          end;
          s.(var)(row,1:length(cstrs)-2) = [cstrs{3:end}];

        end;
      end;
      % DEBUG:      if (mod(nline,100)==0); disp('typecast'); toc; end;
    end;
  end;

  % DEBUG:  disp('end while'); toc;

return;
