1;

 url='http://tds.hycom.org/thredds/dodsC/GOMl0.04/expt_20.1/2004';
 nc = mDataset(url);
 var = 'u';
 dts = cast(getTimes(nc{var}),'double');
 close(nc);

 % Dataset contains daily data - all time stamps should be integral.
 % Therefore, we expect the following check to always return zero:
 dterr = dts - round(dts);
 badix = find( dterr > eps );
 length(badix),
