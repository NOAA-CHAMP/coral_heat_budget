1;

tic,

stns = get_all_station_metadata;

dts = datenum(1987,02,01):(3/24):(datenum(1987,02,02)-eps);

[yrs,mns,dys,hrs] = datevec(dts);

for dix = 1:length(dts)

  dt = dts(dix);

  %http://nomads.ncdc.noaa.gov/thredds/ncss/grid/narr/198701/19870131/narr-a_221_19870131_2100_000.grb?var=Pressure&var=Specific_humidity_height_above_ground&var=Dew_point_temperature&var=Relative_humidity&var=u_wind_height_above_ground&var=v_wind_height_above_ground&var=Downward_longwave_radiation_flux&var=Downward_shortwave_radiation_flux&var=Evaporation&var=Latent_heat_flux&var=Sensible_heat_flux&var=Total_cloud_cover&var=Total_precipitation&var=Upward_long_wave_radiation_flux_surface&var=Upward_short_wave_radiation_flux_surface&latitude=25.01&longitude=-80.38&temporal=all&time_start=1987-01-31T21%3A00%3A00Z&time_end=1987-02-01T00%3A00%3A00Z&time=1987-01-31T21%3A00%3A00Z&vertCoord=0&point=true

  for stix = 1:length(stns.lons)

   if ( stns.codes{stix}(4:5) == 'F1' )

     url = sprintf('http://nomads.ncdc.noaa.gov/thredds/ncss/grid/narr/%04d%02d/%04d%02d%02d/narr-a_221_%04d%02d%02d_%02d00_000.grb?var=Pressure&var=Specific_humidity_height_above_ground&var=Dew_point_temperature&var=Relative_humidity&var=u_wind_height_above_ground&var=v_wind_height_above_ground&var=Downward_longwave_radiation_flux&var=Downward_shortwave_radiation_flux&var=Evaporation&var=Latent_heat_flux&var=Sensible_heat_flux&var=Total_cloud_cover&var=Total_precipitation&var=Upward_long_wave_radiation_flux_surface&var=Upward_short_wave_radiation_flux_surface&latitude=%g&longitude=%g&temporal=point&time=%04d-%02d-%02dT%02d:00:00Z&vertCoord=0&point=true&accept=netcdf', ...
                   yrs(dix), mns(dix), yrs(dix), mns(dix), dys(dix), ...
                   yrs(dix), mns(dix), dys(dix), hrs(dix), ...
                   stns.lats(stix), stns.lons(stix), ...
                   yrs(dix), mns(dix), dys(dix), hrs(dix) );
     %DEBUG:
     url,

     nc = mDataset(url);

     lon(dix,stix) = nc{'longitude'}(1:end,1:end,1:end,1:end);
     lat(dix,stix) = nc{'latitude'}(1:end,1:end,1:end,1:end);
     dsf(dix,stix) = nc{'Downward_shortwave_radiation_flux'}(1:end,1:end,1:end,1:end);
     dlf(dix,stix) = nc{'Downward_longwave_radiation_flux'}(1:end,1:end,1:end,1:end);
     usf(dix,stix) = nc{'Upward_short_wave_radiation_flux_surface'}(1:end,1:end,1:end,1:end);
     ulf(dix,stix) = nc{'Upward_long_wave_radiation_flux_surface'}(1:end,1:end,1:end,1:end);

     close(nc);

   end;

  end;

end;

clear dt dts dx dy dys hrs mns url x yrs;

toc,
