function result = importmetwxtaver(fname)
%function result = importmetwxtaver(fname)
%
% Import an ASCII WXT520 datafile FNAME and return the data struct RESULT
% with the timestamped (hour-averaged) meteorological data from the file.
%
% Last Saved Time-stamp: <Mon 2010-07-26 17:46:01 Eastern Daylight Time gramer>

    result = [];

    if ( ~exist(fname,'file') )
      warning('ImportMetWxtAver:NoFile','Unable to find "%s".',fname);
      return;
    end;
    fid = fopen(fname,'r');
    if ( fid < 0 )
      warning('ImportMetWxtAver:NoFile','Unable to open "%s".',fname);
      return;
    end;
    dat = textscan(fid,'%s %s %f %f %f %f %f %f %f %*f %*f %f %*[^\n]');
    fclose(fid);
    if ( length(dat) < 10 || isempty(dat{1}) )
      warning('ImportMetWxtAver:NoFile','Cannot import "%s".',fname);
      return;
    end;

    alldts = datenum(strcat(dat{1},dat{2}),'yyyymmddHHMMSS');
    [yr,mo,dy,hr,mn,sc] = datevec(alldts);
    dts = unique(datenum(yr,mo,dy,hr,0,0));

    [n,bin] = histc(alldts, [dts ; (floor(dts(end)) + 1)]);

    result.wxt_wspeed.date = dts;
    result.wxt_wspeed.data = mps2kts(grpstats(dat{3},bin));

    result.wxt_wgust.date = dts;
    result.wxt_wgust.data = mps2kts(grpstats(dat{5},bin));


    % Don't forget the vectorial averaging for wind direction!
    [u,v] = spddir_to_uv(dat{3},dat{6});

    result.wxt_u.date = dts;
    result.wxt_u.data = grpstats(u,bin);

    result.wxt_v.date = dts;
    result.wxt_v.data = grpstats(v,bin);

    result.wxt_wdir.date = dts;
    result.wxt_wdir.data = uv_to_dir(result.wxt_u.data,result.wxt_v.data);


    result.wxt_air_t.date = dts;
    result.wxt_air_t.data = grpstats(dat{7},bin);

    result.wxt_relhumid.date = dts;
    result.wxt_relhumid.data = grpstats(dat{8},bin);

    result.wxt_barom.date = dts;
    result.wxt_barom.data = grpstats(dat{9},bin);

    result.wxt_precip.date = dts;
    result.wxt_precip.data = grpstats(dat{10},bin);

return;
