1;

dat=stn.fkeys_hycom_sgs_thermal_diffusivity.data;
dat(dat<0.5 | dat>2000 | ~isfinite(dat)) = [];
dat(dat<15 | dat>2000 | ~isfinite(dat)) = [];

x=dat./(nanmax(dat)+1);
bps=betafit(x),gps=gamfit(x),wps=wblfit(x),

y=betapdf(x,bps(1),bps(2)); figure; maxigraph; plot(dat,y,'.'); titlename('Beta fit');

y=gampdf(x,gps(1),gps(2)); figure; maxigraph; plot(dat,y,'.'); titlename('Gamma fit');

y=wblpdf(x,wps(1),wps(2)); figure; maxigraph; plot(dat,y,'.'); titlename('Weibull fit');
