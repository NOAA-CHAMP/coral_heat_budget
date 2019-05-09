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

% [newq.date,newq.data] = gap_expand(stn.(qtfld).date(goodix),real(stn.(qtfld).data(goodix)));
newq.date = [stn.(qtfld).date(goodix(1)):(1/24):stn.(qtfld).date(goodix(end))]';
newq.data = interp1(stn.(qtfld).date(goodix),real(stn.(qtfld).data(goodix)),newq.date);

%%%% DEBUG??? newq.data = newq.data .* (1/24);

% Truncate initial and final NaN sequences
begix = find(isfinite(newq.data),1);
endix = find(isfinite(newq.data),1,'last');
newq.date = newq.date(begix:endix);
newq.data = newq.data(begix:endix);

% Avoid seasonal bias
begix = find( (get_jday(newq.date)==1 & get_hour(newq.date)==0), 1 );
endix = find( (get_jday(newq.date)==365 & get_hour(newq.date)==23), 1, 'last' );
newq.date = newq.date(begix:endix);
newq.data = newq.data(begix:endix);

% Stupid leap years
badix = find(get_jday(newq.date)==366);
newq.date(badix) = [];
newq.data(badix) = [];


badix = find(ismember(get_year(newq.date),[2005]));
newq.date(badix) = [];
newq.data(badix) = [];


nyrs = length(unique(get_year(newq.date)));

yrhrs = 24*365;
dat = reshape(newq.data,[yrhrs nyrs])';

acc = cumsum(dat,2);
if ( ~doAcc )
  acc = dat - repmat(dat(:,1),[1 yrhrs]);
end;

figure;
maxigraph;
plot(1:(1/24):366-(1/24),nanmedian(acc));
xlim([1 366]);
titlename([stn.station_name ' ' strrep(qtfld,'_','\_')]);
