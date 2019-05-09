1;

doAcc = true;

% qtfld = 'ndbc_sea_t'; doAcc = false;

% qtfld = 'ndbc_hfbulk_heat_flux_term';
% qtfld = 'ndbc_erai_30a_heat_flux_term';
% qtfld = 'ndbc_erai_30a_ww3_fkeys_qe_qtadv';
% qtfld = 'ndbc_erai_30a_ww3_fkeys_qe_dt';
% qtfld = 'benthic_ndbc_erai_30a_ww3_fkeys_qe_dt';
qtfld = 'benthic_ndbc_erai_30a_ww3_fkeys_qe_dt_netqf';

goodix = 1:length(stn.(qtfld).data);
% goodix = ts_boreal_cool(stn.(qtfld));
% goodix = ts_boreal_warm(stn.(qtfld));

[newq.date,newq.data] = gap_expand(stn.(qtfld).date(goodix),stn.(qtfld).data(goodix));


% Truncate initial and final NaN sequences
begix = find(isfinite(newq.data),1);
endix = find(isfinite(newq.data),1,'last');
newq.date = newq.date(begix:endix);
newq.data = newq.data(begix:endix);

% Avoid seasonal bias
begix = find(get_jday(newq.date)==1,1);
endix = find(get_jday(newq.date)==365,1,'last');
newq.date = newq.date(begix:endix);
newq.data = newq.data(begix:endix);

% Avoid (tiny) diurnal bias
begix = find(get_hour(newq.date)==0,1);
endix = find(get_hour(newq.date)==23,1,'last');
newq.date = newq.date(begix:endix);
newq.data = newq.data(begix:endix);

ndys = length(unique(floor(newq.date)));

dat = reshape(newq.data,[24 ndys])';

acc = cumsum(dat,2);
if ( ~doAcc )
  acc = dat - repmat(dat(:,1),[1 24]);
end;

figure;
maxigraph;
plot(0:23,nanmean(acc));
titlename([stn.station_name ' ' strrep(qtfld,'_','\_')]);
