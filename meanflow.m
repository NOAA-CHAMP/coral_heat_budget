function [U, V] = meanflow(LON, LAT, u0, v0)

  if ( ~exist('u0', 'var') || isempty(u0) )
    u0 = 0;
  end;
  if ( ~exist('v0', 'var') || isempty(v0) )
    v0 = 0;
  end;

  U = repmat(u0, size(LON));
  V = repmat(v0, size(LON));

return;
