1;

more off;

if (~exist('smkf1','var'))
    smkf1 = get_station_from_station_name('smkf1');
    smkf1 = load_all_ndbc_data(smkf1);
end;
if (~exist('lonf1','var'))
    lonf1 = get_station_from_station_name('lonf1');
    lonf1 = load_all_ndbc_data(lonf1);
end;

pause;
for yr=[1993:2003,2005:2007];
    plot_spec(smkf1,{'ndbc_wind1_speed','ndbc_air_t','ndbc_sea_t',},[],[],[],[],[],true,[],{@(x)(find(datenum(yr,5,1)<=x.date&x.date<datenum(yr,10,1)))});
    appendtitlename([' ',num2str(yr),' WARM']);
end;

pause;
for yr=[1993:2004,2006:2008];
    plot_spec(smkf1,{'ndbc_wind1_speed','ndbc_air_t','ndbc_sea_t',},[],[],[],[],[],true,[],{@(x)(find(datenum(yr-1,11,1)<=x.date&x.date<datenum(yr,4,1)))});
    appendtitlename([' ',num2str(yr),' COOL']);
end;

pause;
for yr=1993:2011;
    plot_spec(lonf1,{'ndbc_wind1_speed','ndbc_air_t','ndbc_sea_t',},[],[],[],[],[],true,[],{@(x)(find(datenum(yr,5,1)<=x.date&x.date<datenum(yr,10,1)))});
    appendtitlename([' ',num2str(yr),' WARM']);
end;

pause;
for yr=1993:2011;
    plot_spec(lonf1,{'ndbc_wind1_speed','ndbc_air_t','ndbc_sea_t',},[],[],[],[],[],true,[],{@(x)(find(datenum(yr-1,11,1)<=x.date&x.date<datenum(yr,4,1)))});
    appendtitlename([' ',num2str(yr),' COOL']);
end;

clear yr ans

more on
