1;

for callForm = [1 7 9]

doAcc = true;
switch (callForm)
 case 1, qtfld = 'ndbc_sea_t'; doAcc = false;

 case 2, qtfld = 'ndbc_hfbulk_heat_flux_term';
 case 3, qtfld = 'ndbc_erai_30a_heat_flux_term';
 case 4, qtfld = 'benthic_ndbc_erai_30a_net_heat_flux_term';
 case 5, qtfld = 'ndbc_erai_30a_ww3_fkeys_qe_qtadv';
 case 6, qtfld = 'ndbc_erai_30a_ww3_fkeys_qe_dt';
 case 7, qtfld = 'benthic_ndbc_erai_30a_ww3_fkeys_qe_dt';
 case 8, qtfld = 'benthic_ndbc_erai_30a_ww3_fkeys_qe_dt_netqf';

 case 9, qtfld = 'hc_dTdt';
end;

goodix = 1:length(stn.(qtfld).data);

[newq.date,newq.data] = gap_expand(stn.(qtfld).date(goodix),real(stn.(qtfld).data(goodix)));
% newq.date = [stn.(qtfld).date(goodix(1)):(1/24):stn.(qtfld).date(goodix(end))]';
% newq.data = interp1(stn.(qtfld).date(goodix),real(stn.(qtfld).data(goodix)),newq.date);

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


% badix = find(ismember(get_year(newq.date),[2005]));
% newq.date(badix) = [];
% newq.data(badix) = [];


nyrs = length(unique(get_year(newq.date)));

yrhrs = 24*365;
dat = reshape(newq.data,[yrhrs nyrs])';

mdat = nanmean(dat,1);

acc = cumsum(mdat);
%%%% DEBUG ??? 
acc = mdat;
if ( ~doAcc )
  acc = mdat - mdat(1);
  %%%% DEBUG ???   acc = [0 , diff(mdat)];
end;

% figure;
% maxigraph;
% plot(1:(1/24):366-(1/24),acc);
% xlim([1 366]);
% titlename([stn.station_name ' ' strrep(qtfld,'_','\_')]);

% filthr = 72; facc = filter_ts([acc acc],filthr); figure; maxigraph; plot(1:(1/24):366-(1/24),facc(1:8760)); xlim([1 366]); titlename([stn.station_name ' ' strrep(qtfld,'_','\_') ' ' num2str(filthr) 'hlp' ]);

% filthr = 168; facc = sma_ts(acc,filthr); figure; maxigraph; plot(1:(1/24):366-(1/24),facc(1:8760)); xlim([1 366]); titlename([stn.station_name ' ' strrep(qtfld,'_','\_') ' ' num2str(filthr) 'hsma' ]);

ann = grpstats(newq.data,(24*get_yearday(newq.date)));
if ( ~doAcc )
  ann = ann - ann(1);
else
  ann = cumsum(ann);
end;

figure;
maxigraph;
plot(1:(1/24):366-(1/24),ann);
xlim([1 366]);
titlename([stn.station_name ' ' strrep(qtfld,'_','\_')]);

end;
