1;
%% SCRIPT to estimate the error magnitude for gamma (absorbed net insolation)
%
% Last Saved Time-stamp: <Thu 2012-02-23 14:51:42  lew.gramer>


%% Attenuation coefficient for combined PAR/UV light
% Kd = 0.10;
%% Fraction of sea floor at site which is highly reflective
% sandfrac = station_sand_cover(stn);


dts = [datenum(2012,1,1):(1/24):datenum(2012,12,32)]';

if ( exist('stn','var') )
  lon = stn.lon;
  lat = stn.lat;

  hfld = ['tpxo_i_depth'];
  mhfld = ['mean_' hfld];
  if ( isfield(stn,mhfld) )
    h = nanmean(stn.(mhfld).data);
  elseif ( isfield(stn,'depth') )
    h = stn.depth;
  end;
end;
if ( ~exist('lon','var') )
  lon = -80.38;
  lat = +25.01;
end;
if ( ~exist('h','var') )
  h = 2;
end;

Ppen = 0.501 + 0.048;

for Kd=[0.05:0.15:0.80];
  % for sandfrac = [0.1:0.1:0.4];
  % LONF1 sandfrac=0.50, MLRF1 sandfrac=0.30
  for sandfrac = [0.3:0.1:0.5];
    % Sea-bottom mean reflectance
    Ab = (0.40 * sandfrac) + (0.07 * (1-sandfrac));
    [yds,yrs] = get_yearday(dts);
    [theta,ig] = soradna1(yds,yrs,-lon,lat);
    sun_angle_correction = secd(90-theta);
    tau = exp(-Kd .* h .* sun_angle_correction);
    tau_scatter = exp(-Kd .* h .* secd(90-45));
    % gam = ( 1 - (Ab.*tau.*(1 - tau)) );
    gam = (1 - Ppen) + (Ppen.*(1 - tau)) + (Ppen.*tau.*Ab.*(1-tau_scatter));

    gam(theta < 5) = 0;
    %???DEBUG:    gam(theta < 20) = 0;
    %???DEBUG:
    gam(gam<0.01) = [];

    % disp(['Kd,Sand: ', num2str([Kd,sandfrac])]);
    % nansummary(gam);

    % Standard deviation of "error"
    ersd = (max(gam)-min(gam))/6;
    disp(['Kd,Sand: ',num2str([Kd,sandfrac]),'...',num2str(prctile(gam,7),'%.2f'),'-',num2str(prctile(gam,93),'%.2f'),'...',num2str(ersd,'%.2f'),'=',num2str(nanmean(100*ersd./gam),'%.0f'),'% of gamma']);
  end;
end;

