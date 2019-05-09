1;

sfld = 'ndbc_sea_t';
hcdTdt = 'ndbc_erai_ww3_30a_ww3_avhrr_hc_dTdt';

cmpfld = hcdTdt;
errfld = [hcdTdt,'_err'];
[stn.t.data,stn.t.date] = grp_ts(stn.(sfld).data,stn.(sfld).date,@floor,@nanmean);
[stn.q.data,stn.q.date] = grp_ts(stn.(cmpfld).data,stn.(cmpfld).date,@floor,@sum,24);
[stn.e.data,stn.e.date] = grp_ts(stn.(errfld).data,stn.(errfld).date,@floor,@sum,24);

stn.model_t = ts_op(stn.t,stn.q,'+');
stn.model_t_err_lo = ts_op(stn.model_t,stn.e,'-');
stn.model_t_err_hi = ts_op(stn.model_t,stn.e,'+');

stn.model_err = ts_op(stn.t,stn.model_t,'-');
fmg;
plot_ts(stn.model_err);
titlename(strrep([upper(stn.station_name),' Error: ',sfld,' daily mean vs. ',cmpfld],'_','\_'));

stn.model_invalid = ts_op(stn.model_err,stn.e,@(x,y)(abs(x)>abs(y)));
stn.model_invalid = subset_ts(stn.model_invalid,find(stn.model_invalid.data));
[ig,stn.model_t_err_lo_invalid] = intersect_tses(stn.model_invalid,stn.model_t_err_lo);
[ig,stn.model_t_err_hi_invalid] = intersect_tses(stn.model_invalid,stn.model_t_err_hi);

fmg;
% plot_ts(stn.t,'k.-','LineWidth',3,stn.model_t,'r.-');
% plot_ts(stn.model_t_err_lo_invalid,'r^',stn.model_t_err_hi_invalid,'rv');
plot_ts(stn.t,'k.','LineWidth',3);
plot_ts(stn.model_t_err_lo,'b-',stn.model_t_err_hi,'r-');
legend('T_s','\Sigma_2_4_h[\partial_tT_s \pm\sigma]');
titlename(strrep([upper(stn.station_name),' ',sfld,' daily mean vs. ',cmpfld],'_','\_'));

fmg; hist(stn.model_err.data,100);
