function [U, V] = eddy(LON, LAT, U, V, maxVr, vprof, lon, lat, radiuskm)
%function [U, V] = eddy(LON, LAT, U, V, maxVr, vprof, lon, lat, radiuskm)
% 
% Generate a test eddy on the given LON/LAT grid. Generates radial velocities
% Vr up to a maximum maxVr [m/s], centered at 'lon' and 'lat', and perturbs
% the mean field [U,V] (optionally empty matrices) with those Vrs. Profile of
% Vr as function of radial distance is determined by 'vprof' string argument:
%     'sech' - (DEFAULT) hyperbolic secant profile centered at middle radius
%     'sine' - sine-bell profile centered at middle radius 'radiuskm/2'
%     'solid' - solid-body rotation with maxVr at the outer edge 'radiuskm'
% 
% If 'maxVr' < 0, a COUNTER-CLOCKWISE eddy (cyclonic in NHem) is generated.
% 

  radiusm = radiuskm * 1e3;

  if ( isempty(U) )
    U = repmat(nan, size(LON));
  end;
  if ( isempty(V) )
    V = repmat(nan, size(LON));
  end;

  sech_factor = 1;
  if ( isempty(vprof) )
    vprof = 'sech';
  elseif ( isnumeric(vprof) )
    sech_factor = vprof(1);
    vprof = 'sech';
  end;

  % Slow-but-readable loop over all nearby grid points (for now)
  % for ix = 1:numel(U(boxix))
  radiuslat = (radiuskm*1.1) / 111.31;
  radiuslon = radiuslat * cosd(lat);
  boxix = find(abs(LON - lon) < radiuslon & abs(LAT - lat) < radiuslat);
  for idx = 1:length(boxix)

    ix = boxix(idx);
    [rngkm, brgE] = sw_dist([lat LAT(ix)], [lon LON(ix)], 'km');
    % Defaults: nautical miles, and degrees CCW from East. WTF, CSIRO??
    rng = rngkm * 1e3;
    if (brgE <= 90); brg = 90 - brgE; else; brg = 450 - brgE; end;

    % Innermost 1km of the eddy is quiescent?
    if ( 1e3 <= rng && rng <= radiusm )
      switch(vprof)
       case { 'solid' },
        % Solid-body rotation
        Vr = maxVr*(rng/radiusm);
       case { 'sin', 'sine' },
        % Maximum velocity at 1/2 radius - sine-bell profile
        Vr = maxVr*sin(pi*(rng/radiusm));
       case { 'sech' },
        % Maximum velocity at 1/2 radius - "sharp" sech profile
        x = sech_factor * 10 * ( (rng/radiusm) - 0.5 );
        Vr = maxVr*sech(x);
       otherwise,
        error('Eddy:UnknownProfile', ...
              'Unknown v-profile type-string "%s"!', vprof);
      end;

      U(ix) = U(ix) + cosd(brg)*Vr;
      V(ix) = V(ix) - sind(brg)*Vr;

    end;

  end;

return;
