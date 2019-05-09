1;
%%%% SCRIPT to create harmonic tidal analysis for JMR MS. from Ch. 2

s = smkf1;
%s = lonf1;

disp(s.station_name);

[tin,xin] = gap_expand(s.ndbc_tide_m.date,s.ndbc_tide_m.data);

args = {'start time',tin(1),'latitude',s.lat};
[tres,xout] = t_tide(xin,args{:},'output',fullfile(get_thesis_path('../src'),[mfilename,'-',s.station_name,'.log']));
%[tres,xout] = t_tide(xin,args{:},'output','none');
tres.headers = char({'amp','amp_err','pha','pha_err'});

[ig,sortix] = sort(tres.tidecon(:,1),1,'descend');
tres.namesorted = tres.name(sortix,:);
tres.freqsorted = tres.freq(sortix,:);
tres.tidesorted = tres.tidecon(sortix,:);
clear ig sortix

%maxix=find(tres.tidesorted(:,1)>prctile(tres.tidesorted(:,1),90));
% Tide reported to 0.01 feet == 0.003 m
maxix=find(tres.tidesorted(:,1)>0.006);
for ix=maxix(:)';
  if ( tres.freqsorted(ix) > 1/36 )
    freqstr = sprintf('%6.2f h',1/tres.freqsorted(ix));
  else
    freqstr = sprintf('%6.2f d',1/tres.freqsorted(ix)/24);
  end;
  disp(sprintf('%5s %s  %f  %f  %5.1f  %4.1f',tres.namesorted(ix,:),freqstr,tres.tidesorted(ix,:)));
end;
clear ix freqstr

r.date=tin; r.data=xin; m.date=tin; m.data=xout;
[B,Stats,fh] = scatter_fit_ts(m,r,[],[],'Real','T_TIDE',[],'resid');
figure(fh(2)); hold on;
plot(s.ndbc_wind1_speed.date,s.ndbc_wind1_speed.data./100,'r-');
clear r m

resid.date = Stats.datenums;
resid.data = Stats.resid;

scatter_fit_ts(s.ndbc_wind1_speed,resid,[],[],'Wind','h_t_i_d_e Resid');
scatter_fit_ts(s.ndbc_wind1_speed,resid,@(x)(find(x.data>30)),[],'Wind >30kts','h_t_i_d_e Resid');
scatter_fit_ts(s.ndbc_barom,resid,[],[],'Barom','h_t_i_d_e Resid');

x.resid = resid;
plot_spec(x,'resid',[],[],[],[],[],true,true);
x=[]; clear x;

s=[]; clear s;
