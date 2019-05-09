function [U, V] = randflow(LON, LAT, U, V, u0)
%function [U, V] = randflow(LON, LAT, U, V, u0)
%
% Add stochastic noise to a flow field: u0 is the desired magnitude of the
% standard deviation for the normally distributed noise.
%

  if ( ~exist('u0', 'var') || isempty(u0) )
    u0 = 0.1 * std(sqrt((U.^2) + (V.^2)));
  end;

  U = U + ( 2*u0*(randn(size(U)) - 0.5) );
  V = V + ( 2*u0*(randn(size(V)) - 0.5) );

return;
