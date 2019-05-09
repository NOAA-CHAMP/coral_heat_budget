1;

if ( ~exist('stn','var') || isempty(stn) )
    stn = get_station_from_station_name('mlrf1');
end;

if ( ~isfield(stn,'ndbc_wind1_speed') )
    stn = load_all_ndbc_data(stn);
end;

% hrs = [15:19];
hrs = 19;
mos = 1:2:12;
dy = 15;
Ws = 0:30;
Cs = 0.00:0.20:1.00;

cs = {'ko','kx','bo','bx','ro','rx','mo','mx','co','cx','go','gx'};

for hrix=1:numel(hrs)
    hr = hrs(hrix);

    fmg;
    for moix=1:numel(mos)
        mo = mos(moix);
        dts = datenum(2011,mo,dy,hr,0,[0:numel(Ws)-1]);
        stn.W.date = dts';
        stn.W.data = Ws';

        subplot_tight(ceil(numel(mos)/2),2,moix);
        hold on;
        for Cix=1:numel(Cs)
            C = Cs(Cix);
            stn.C = stn.W; stn.C.data(:) = C;
            stn = station_bulk_albedo(stn,'alb','W','C');
            plot(stn.W.data,stn.alb.data,cs{Cix});
        end;
        xlabel('Wind speed [kts]');
        ylabel(['\alpha_S_W']);
        ylim([0.0,0.2]);
        legh = legend(strcat('C=',num2str(Cs')), 'Location','NorthEast');
        xlabel(legh,'Cloud fraction');
        title(datestr(dts(1)));
    end;

    suptitlename(['Bulk Albedo Sensitivity Test: Day-Hour ',num2str(hr)]);
end;

clear hr* mo* dy* W* C* cs dt* ans
