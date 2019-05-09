1;

if ( ~exist('stn','var') )
  stn = get_station_from_station_name('mlrf1');
  stn = load_all_ndbc_data(stn);
  stn = load_station_data(stn);
  stn = get_satpar_insol(stn);
  stn = get_erai_station(stn);

  stn = load_cleaning_dates(stn);

  stn = station_insol_to_par(stn,'erai_dsrf','erai_actual_par','erai_actual_parW');
  stn = station_insol_to_par(stn,'sat_insol_in','sat_actual_par','sat_actual_parW');
  stn = station_par_to_insol(stn,'bic_surf_par','bic_surf_dsrf','bic_surf_usrf','bic_surf_srf');

  [ins,sat,era] = intersect_tses(stn.bic_surf_dsrf,stn.sat_insol_in,stn.erai_dsrf);
end;

yrs = sprintf(' (%s to %s)',datestr(ins.date(1)),datestr(ins.date(end)));

fmg;
climins = grpplot_ts(ins,@get_jday,@nanmean,0,'k');
climsat = grpplot_ts(sat,@get_jday,@nanmean,0,'b');
climera = grpplot_ts(era,@get_jday,@nanmean,0,'r');
legend('In situ calculated','TOMS GSIP','ERA-Interim', 'Location','South');
titlename([upper(stn.station_name),' Q_S_W^I climatology',yrs]);

eramins = ts_op(climera,climins,'-');
satmera = ts_op(climsat,climera,'-');
satmins = ts_op(climsat,climins,'-');

fmg;
plot_ts(eramins,satmins,satmera);
legend('ERA-Interim - In situ calc.','TOMS GSIP - In situ calc.','TOMS GSIP - ERA-Interim', 'Location','South');
titlename([upper(stn.station_name),' Q_S_W^I climatologies',yrs]);


fmg;
dlyins = grpplot_ts(ins,@floor,@nanmean,0,'k');
dlysat = grpplot_ts(sat,@floor,@nanmean,0,'b');
dlyera = grpplot_ts(era,@floor,@nanmean,0,'r');
datetick3('x',2,'keeplimits');
legend('In situ calculated','TOMS GSIP','ERA-Interim', 'Location','South');
titlename([upper(stn.station_name),' Q_S_W^I daily mean',yrs]);

yl = ylim(gca); dy=(yl(2)-yl(1)); yl(1)=yl(2)-(dy/10);
for ix=1:length(stn.cleaning_date); arrow([stn.cleaning_date(ix),yl(2)],[stn.cleaning_date(ix),yl(1)]); end;


dlyeramins = ts_op(dlyins,dlyera,'-');
dlysatmera = ts_op(dlysat,dlyera,'-');
dlysatmins = ts_op(dlysat,dlyins,'-');
fmg;
plot_ts(dlyeramins,dlysatmins,dlysatmera);
legend('ERA-Interim - In situ calc.','TOMS GSIP - In situ calc.','TOMS GSIP - ERA-Interim', 'Location','South');
titlename([upper(stn.station_name),' Q_S_W^I',yrs]);

yl = ylim(gca); dy=(yl(2)-yl(1)); yl(1)=yl(2)-(dy/10);
for ix=1:length(stn.cleaning_date); arrow([stn.cleaning_date(ix),yl(2)],[stn.cleaning_date(ix),yl(1)]); end;


[B,Stats]=scatter_fit_ts(dlyera,dlyins,[],[],'ERA-Interim Q_S_W^I','In situ calc. Q_S_W^I',[],[],true),
[B,Stats]=scatter_fit_ts(dlysat,dlyins,[],[],'TOMS GSIP Q_S_W^I','In situ calc. Q_S_W^I',[],[],true),
