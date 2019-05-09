1;

PFX = '';

euleru = [PFX 'quasi_eulerian_u'];
eulerv = [PFX 'quasi_eulerian_v'];

UdotdelT = [PFX 'advected_heat'];

if ( isfield(stn,UdotdelT) )
  stn = rmfield(stn,UdotdelT);
end;

stn = station_advect_field(stn,UdotdelT,...
                           euleru,eulerv,...
                           'avhrr_weekly_sst');
