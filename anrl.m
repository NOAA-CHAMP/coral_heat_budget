1;

% Analyze heat flux run-lengths - numbers of hours between zero crossings
% and full forcing periods, i.e., changes from positive to negative, then
% back to positive surface heat flux. Assumes STN contains a heat budget.

clear ans st xix upix badix dt fs

% % Surface heat flux: Q0 = QSW+QLW+QLH+QSW+QRH
% st.q=stn.simple_ndbc_erai_erai_30a_net_flux;
% Net vertical heat flux (incl. benthic flux, lost insolation): Q0+Qb
st.q=stn.b_ndbc_erai_erai_30a_net_flux;

st = verify_variable(st,'q_24_h_avg');

% Sign (ocean cooling or heating) of surface heat flux
st.s.date=st.q.date;
st.s.data=sign(st.q.data);
%st=filter_gaps(st,'q','s',(3.1/24),0);

% Time series of sign differences in ST.Q
st.ds.date=st.s.date(2:end);
st.ds.data=diff(st.s.data);
%st=filter_gaps(st,'q','ds',(3.1/24),0);

% Time series of run-lengths between zero-crossings in ST.Q
xix=[1;find(abs(st.ds.data)>eps)+1;numel(st.q.date)];
st.rl.date=st.q.date(xix(2:end));
st.rl.data=diff(st.q.date(xix))*24;
%st=filter_gaps(st,'q','rl',(3.1/24),0);
st=filter_gaps(st,'q','rl',(3.1/24),1);
% Quality control: No cold snap lasts longer than 7 days
badix=find(st.rl.data>(24*7));
if ( ~isempty(badix) )
  numel(badix), st.rl.data(badix([1,end])), datestr(st.rl.date(badix([1,end]))),
  st.rl.date(badix)=[];
  st.rl.data(badix)=[];
end;

% Time series of run-lengths of warming and cooling periods in ST.Q, resp.
six=find(ismember(st.s.date,st.rl.date));
st.prl.date=st.rl.date(st.s.data(six-1)>eps);
st.prl.data=st.rl.data(st.s.data(six-1)>eps);
st.nrl.date=st.rl.date(st.s.data(six-1)<eps);
st.nrl.data=st.rl.data(st.s.data(six-1)<eps);

% Time series of run-lengths between UPWARD zero-crossings in ST.Q
upix=[1;find(st.ds.data>eps)+1;numel(st.q.date)];
st.url.date=st.q.date(upix(2:end));
st.url.data=diff(st.q.date(upix))*24;
%st=filter_gaps(st,'q','url',(3.1/24),0);
st=filter_gaps(st,'q','url',(3.1/24),1);
% First forcing cycle is almost always partial
st.url.date(1)=[]; st.url.data(1)=[]; 
% Quality control
badix=find(st.url.data>(24*7));
if ( ~isempty(badix) )
  numel(badix), st.url.data(badix([1,end])), datestr(st.url.date(badix([1,end]))),
  st.url.date(badix)=[];
  st.url.data(badix)=[];
end;

% % fmg; plot(st.rl.date,st.rl.data); datetick3;
% % fmg; boxplot_ts(st.rl,[],'mean',true);
% fmg; plot(st.url.date,st.url.data); datetick3;
% fmg; boxplot_ts(st.url,[],'mean',true);

% Dates of some cold snaps in the Keys
% dt=datenum(1989,03,08); f=@(x)(find(dt-14<=x.date&x.date<=dt+14));
% dt=datenum(1989,10,13); f=@(x)(find(dt-30<=x.date&x.date<=dt+30));
% dt=datenum(1991,12,20); f=@(x)(find(dt-14<=x.date&x.date<=dt+14));
% dt=datenum(1993,11,01); f=@(x)(find(dt-14<=x.date&x.date<=dt+14));
% dt=datenum(1993,12,25); f=@(x)(find(dt-14<=x.date&x.date<=dt+14));
% dt=datenum(2000,01,31); f=@(x)(find(dt-14<=x.date&x.date<=dt+14));
dt=datenum(2000,02,11); f=@(x)(find(dt-45<=x.date&x.date<=dt+45));
% dt=datenum(2007,12,31); f=@(x)(find(dt-14<=x.date&x.date<=dt+14));

fmg; plot_ts(subset_ts(st.q,f),'b.-',subset_ts(st.q_24_h_avg,f),'b-'); titlename('Flux');
plot_ts(subset_ts(subset_ts(st.q,xix),f),'ko','MarkerSize',8);
% fmg; plot_ts(subset_ts(st.s,f)); titlename('Flux sign');
% fmg; plot_ts(subset_ts(st.ds,f)); titlename('Flux sign change');
fmg; plot_ts(subset_ts(st.rl,f),'b.-'); titlename('Zero-crossing runlengths');
fmg; plot_ts(subset_ts(st.nrl,f),'b.-',subset_ts(st.prl,f),'r.-'); titlename('Warming/cooling period runlengths');
fmg; plot_ts(subset_ts(st.url,f),'b.-'); titlename('Warming zero-crossing runlengths');
