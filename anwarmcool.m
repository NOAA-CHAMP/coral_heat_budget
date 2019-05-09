1;

%fld='ndbc_erai_erai_30a_net_flux';
fld='hourlymisst_erai_erai_erai_30a_net_flux';
trmfld = [fld,'_term'];

[q0,qt] = intersect_tses(stn.(fld),stn.(trmfld));

gapix=unique([1 ; ...
    find(diff(q0.date)>(3.01/24))+1 ; ...
    length(q0.date)+1]);

perlen=[];
perix=[];

for gapixix=2:numel(gapix)
    ixix=gapix(gapixix-1):(gapix(gapixix)-1);
    ix=[find(diff(sign(q0.data(ixix)))~=0)+1];
    ix=unique(ixix([1;ix]));
    perlen=[ perlen(:) ; diff(q0.date(ix)) ];
    perix=[ perix(:) ; ix(1:end-1)' ];
end;

percum=repmat(0,size(perix));
for perixix=2:numel(perix)
    percum(perixix-1) = sum(qt.data(perix(perixix-1):perix(perixix)-1));
end;
percum(end) = sum(qt.data(perix(perixix):end));

fmg; hist(percum,100); titlename('Distribution of cumulative continuous cooling/warming'); xlabel('\Sigma_c_/_wQ_0/\rhoC_ph [K]');
fmg; hist(perlen,100); titlename('Distr. period lengths of continuous cooling/warming'); xlabel('Days');
fmg; hist(perlen(percum>0),100); titlename('Distr. period lengths of continuous warming'); xlabel('Days');
fmg; hist(perlen(percum<0),100); titlename('Distr. period lengths of continuous cooling'); xlabel('Days');

stn=verify_variable(stn,[fld,'_term_1_d_sum']);

fh=fmg;
plot_ts(stn.([fld,'_term_1_d_sum']),stn.ndbc_air_t,stn.ndbc_sea_t);
legend('\Sigma_1_dQ_0(\gamma)','T_a','T_s', 'Location','South');
for ix=perix(perlen>2)'
    if ( ~ishandle(fh) )
        break;
    end;
    figure(fh);
    xlim(q0.date(ix)+[-5,+5]);
    datetick3;
    titlename([stn.station_name,'.',strrep(fld,'_','\_'),' ',datestr(q0.date(ix))]);
    pause;
end;
