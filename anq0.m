1;

error('Run-once testing script!');

% stnm = 'smkf1';
stnm = 'mlrf1';

station = load_station_data(stnm);

% station = qa_ts(station, 'sea_t', [], true);

dat = read_all_insolation_csvs(stnm);
station.sat_par.date = dat.date;
station.sat_par.data = dat.par;
station.sat_par_watt.date = dat.date;
station.sat_par_watt.data = dat.parW;
station.sat_insol_in.date = dat.date;
station.sat_insol_in.data = dat.sd;
station.sat_insol_out.date = dat.date;
station.sat_insol_out.data = dat.sd;
station.sat_net_insol.date = dat.date;
station.sat_net_insol.data = dat.sd - dat.su;
station.sat_longwave_in.date = dat.date;
station.sat_longwave_in.data = dat.ld;
station.sat_longwave_out.date = dat.date;
station.sat_longwave_out.data = dat.ld;
station.sat_net_longwave.date = dat.date;
station.sat_net_longwave.data = dat.ld - dat.lu;
dat = []; clear dat;

if ( ~strcmpi(stnm, 'smkf1') )
  smkf1 = load_station_data('smkf1');
  station.air_t_dewp = smkf1.air_t_dewp;
  smkf1 = []; clear smkf1;
end;

station = verify_variable(station, 'wind1_speed_1_day_average');
station = verify_variable(station, 'air_t_1_day_average');
station = verify_variable(station, 'air_t_dewp_1_day_average');
station = verify_variable(station, 'barom_surf_1_day_average');
station = verify_variable(station, 'sea_t_1_day_average');
station = verify_variable(station, 'sea_t_5_hour_findiff');
station = verify_variable(station, 'licor_surf_par_1_day_average');
station = verify_variable(station, 'licor_surf_par_1_day_lowpass');
station = verify_variable(station, 'licor_surf_par_1_day_highpass');
station = verify_variable(station, 'sat_par_1_day_average');
station = verify_variable(station, 'sat_par_watt_1_day_average');
station = verify_variable(station, 'sat_insol_in_1_day_average');
station = verify_variable(station, 'sat_insol_in_1_day_lowpass');
station = verify_variable(station, 'sat_insol_in_1_day_highpass');
station = verify_variable(station, 'sat_insol_out_1_day_average');
station = verify_variable(station, 'sat_net_insol_1_day_average');
station = verify_variable(station, 'sat_longwave_in_1_day_average');
station = verify_variable(station, 'sat_longwave_out_1_day_average');
station = verify_variable(station, 'sat_net_longwave_1_day_average');


% multiplot_station(station, { 'licor_surf_par','licor_surf_par_1_day_average', ...
%                     'sat_insol_in','sat_insol_in_1_day_average' }, ...
%                   'Light data - Molasses Reef', 'Date', ...
%                   {'PAR_i', '\mu_1_d(PAR_i)', 'Insol_s', '\mu_1_d(Insol_s)'});

hiix1 = find(station.licor_surf_par.data > 200);
hiix2 = find(station.sat_insol_in.data > 200);
[ix1,ix2] = intersect_dates(station.licor_surf_par.date(hiix1), ...
                            station.sat_insol_in.date(hiix2));
[B,Stats,fh] = scatter_fit(station.licor_surf_par.data(hiix1(ix1)), ...
                           station.sat_insol_in.data(hiix2(ix2)), ...
                           'PAR_i', 'Insol_s');

[ix1,ix2] = intersect_dates(station.sat_par_1_day_average.date, ...
                            station.sat_insol_in_1_day_average.date);
[B,Stats,fh] = scatter_fit(station.sat_par_1_day_average.data(ix1(1:24:end)), ...
                           station.sat_insol_in_1_day_average.data(ix2(1:24:end)), ...
                           'PAR_s', 'Insol_s');

[ix1,ix2] = intersect_dates(station.sat_par_watt_1_day_average.date, ...
                            station.sat_insol_in_1_day_average.date);
[B,Stats,fh] = scatter_fit(station.sat_par_watt_1_day_average.data(ix1(1:24:end)), ...
                           station.sat_insol_in_1_day_average.data(ix2(1:24:end)), ...
                           'PAR_s_W', 'Insol_s');

station = station_heat_flux(station,'wind1_speed','air_t','air_t_dewp',...
                            'barom','sea_t','sat_net_insol','sat_net_longwave');
station = verify_variable(station, 'relhumid_1_day_average');

% multiplot_station(station, {'sea_t','air_t','sea_t_5_hour_findiff','relhumid','heat_flux_term','latent_heat_flux','sensible_heat_flux',});
multiplot_station(station, {'sea_t','air_t','relhumid','heat_flux_term','latent_heat_flux','sensible_heat_flux',});




[ixiw,ixwi] = intersect_dates(station.sat_longwave_in_1_day_average.date, ...
                              station.wind1_speed_1_day_average.date);
[ixis,ixsi] = intersect_dates(station.sat_longwave_in_1_day_average.date, ...
                              station.licor_surf_par_1_day_average.date);
[ixia,ixai] = intersect_dates(station.sat_longwave_in_1_day_average.date, ...
                              station.air_t_1_day_average.date);
[ixir,ixri] = intersect_dates(station.sat_longwave_in_1_day_average.date, ...
                              station.air_t_dewp_1_day_average.date);
[ixsa,ixas] = intersect_dates(station.licor_surf_par_1_day_average.date, ...
                              station.air_t_1_day_average.date);
[ixaw,ixwa] = intersect_dates(station.air_t_1_day_average.date, ...
                              station.wind1_speed_1_day_average.date);
[ixar,ixra] = intersect_dates(station.air_t_1_day_average.date, ...
                              station.air_t_dewp_1_day_average.date);

[ixrw,ixwr] = intersect_dates(station.air_t_dewp_1_day_average.date, ...
                              station.wind1_speed_1_day_average.date);
[ixrs,ixsr] = intersect_dates(station.air_t_dewp_1_day_average.date, ...
                              station.licor_surf_par_1_day_average.date);

ixi = intersect(ixis,intersect(ixia,ixir));
ixw = intersect(ixwi,intersect(ixwa,ixwr));
ixa = intersect(ixai,intersect(ixas,ixar));
ixs = intersect(ixsi,intersect(ixsa,ixsr));
ixr = intersect(intersect(ixri,ixrs),intersect(ixra,ixrw));

regresp = station.sat_longwave_in_1_day_average.data(ixi(1:24:end))';
regdata = [ ...
    station.air_t_1_day_average.data(ixa(1:24:end)) ; ...
    station.air_t_dewp_1_day_average.data(ixr(1:24:end)) ; ...
    station.wind1_speed_1_day_average.data(ixw(1:24:end)) ; ...
    station.licor_surf_par_1_day_average.data(ixs(1:24:end)) ; ...
          ]';
% stats = regstats(regresp, regdata, 'quadratic');
stats = regstats(regresp, regdata, 'linear');
% Save a little memory - remove projection matrix
stats.hatmat = []; stats = rmfield(stats, 'hatmat');


y = station.air_t_dewp_1_day_average.data(ixr(1:24:end))';
x = station.sat_longwave_in_1_day_average.data(ixi(1:24:end))';
[Bli,Statsli,hli] = scatter_fit(x,y,'L_i','Tdew');

% Clean up workspace a bit
clear ixiw ixwi ixis ixsi ixia ixai ixsa ixas ixaw ixwa ixir ixar ixri ixra ixi ixw ixa ixs ixr regresp regdata;
