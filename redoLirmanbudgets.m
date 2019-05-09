1;

for cst={'bnpin','bnpmi','bnpon','bnpnn','bnppa'};
    % for st={'seatemp','hourly_misst'};
    for st={'hourly_misst'};
        stn=optimize_station_heat_budget(cst{:},'erai','avhrr_weekly','erai','tpxo_tide','erai',st{:},[],'seatemp');
        axis([datenum(2010,1,[1,28]),6,30]); datetick3;
        stn=[]; clear stn; pack;
    end;
end;
clear cst st;
