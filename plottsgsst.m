1;

if ( ~exist('mlrf1','var') )
    mlrf1 = load_all_ndbc_data([],'mlrf1');
    mlrf1 = get_avhrr_weekly_field(mlrf1);
    mlrf1 = station_optimal_isobath_orientation(mlrf1);
    mlrf1 = get_ngdc_bathy_station(mlrf1);
end;
if ( ~exist('tsg','var') )
    tsg = antsg;
end;

if ( ~exist('ds','var') )
    ds = distance_wgs84(mlrf1.lat,mlrf1.lon,tsg.lat,tsg.lon);
end;

jd1 = datenum(0,7,15);
jd2 = datenum(0,9,15);

% Side-length (in degrees) of square of interest around station site
dc = 0.12;

%ix=find(ismember(get_jday(tsg.date),[jd1:jd2]) & ds<(dc*111));
ix=find(ismember(get_jday(tsg.date),[jd1:jd2]) & ds<(dc*111) & abs(mod(mlrf1.isobath_orientation+90,180)-mod(tsg.hdg,180))<45);
fmg;
% [cs,ch]=contour(mlrf1.ngdc_92m_bathy.lon,mlrf1.ngdc_92m_bathy.lat,mlrf1.ngdc_92m_bathy.field,[0:-2:-10,-20,-30]);
% clabel(cs,ch,[-2 -6 -10 -20 -30]);
% axis([mlrf1.lon-dc,mlrf1.lon+dc,mlrf1.lat-dc,mlrf1.lat+dc]); axis square;
% axes; hold on; set(gca,'Color','none');
scatter(tsg.lon(ix),tsg.lat(ix),[],get_year(tsg.date(ix)));
set_pcolor_cursor; colorbar('East'); plot(mlrf1.lon,mlrf1.lat,'kp'); linkaxes;
axis([mlrf1.lon-dc,mlrf1.lon+dc,mlrf1.lat-dc,mlrf1.lat+dc]); axis square;
titlename('TSG Datapoints by Year');

%for yr=[1999,2000,2004,2006]
for yr=1995:2011
    yrix=find(get_year(tsg.date(ix)) == yr);
    if ( ~isempty(yrix) )
        yrix=ix(yrix);
        fmg;
        % [cs,ch]=contour(mlrf1.ngdc_92m_bathy.lon,mlrf1.ngdc_92m_bathy.lat,mlrf1.ngdc_92m_bathy.field,[0:-2:-10,-20,-30]);
        % clabel(cs,ch,[-2 -6 -10 -20 -30]);
        % axis([mlrf1.lon-dc,mlrf1.lon+dc,mlrf1.lat-dc,mlrf1.lat+dc]); axis square;
        % axes; hold on; set(gca,'Color','none');
        % % scatter(tsg.lon(yrix),tsg.lat(yrix),[],tsg.sst(yrix));
        mindsst = -1; maxdsst = +1;
        scatter(tsg.lon(yrix),tsg.lat(yrix),[],sign(cosd(stn.isobath_orientation+90-tsg.hdg(yrix))).*tsg.dsst(yrix));
        set_pcolor_cursor; colorbar('East'); plot(mlrf1.lon,mlrf1.lat,'kp'); linkaxes;
        axis([mlrf1.lon-dc,mlrf1.lon+dc,mlrf1.lat-dc,mlrf1.lat+dc,mindsst,maxdsst,mindsst,maxdsst]); axis square;
        titlename(['TSG SST Year ',num2str(yr)]);
    end;
end;


for yr=[2001 2006]
    yrix=find(get_year(tsg.date(ix)) == yr);
    if ( ~isempty(yrix) )
        yrix=ix(yrix);
        fmg;
        % [cs,ch]=contour(mlrf1.ngdc_92m_bathy.lon,mlrf1.ngdc_92m_bathy.lat,mlrf1.ngdc_92m_bathy.field,[0:-2:-10,-20,-30]);
        % clabel(cs,ch,[-2 -6 -10 -20 -30]);
        % axis([mlrf1.lon-dc,mlrf1.lon+dc,mlrf1.lat-dc,mlrf1.lat+dc]); axis square;
        % axes; hold on; set(gca,'Color','none');
        mindsst = 22; maxdsst = 32;
        scatter(tsg.lon(yrix),tsg.lat(yrix),[],tsg.sst(yrix));
        % mindsst = -1; maxdsst = +1;
        % scatter(tsg.lon(yrix),tsg.lat(yrix),[],sign(cosd(stn.isobath_orientation+90-tsg.hdg(yrix))).*tsg.dsst(yrix));
        set_pcolor_cursor; colorbar('East'); plot(mlrf1.lon,mlrf1.lat,'kp'); linkaxes;
        axis([mlrf1.lon-dc,mlrf1.lon+dc,mlrf1.lat-dc,mlrf1.lat+dc,mindsst,maxdsst,mindsst,maxdsst]); axis square;
        titlename(['TSG SST Year ',num2str(yr)]);
    end;
end;
