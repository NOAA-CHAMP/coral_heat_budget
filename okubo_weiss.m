function [CRLR, DIVR, STRN, STRS, DFRM, OWP] = okubo_weiss(U, V, dx, dy)
%function [CRLR, DIVR, STRN, STRS, DFRM, OWP] = okubo_weiss(U, V, dx, dy)
%
% Calculate curl, div, normal and shear strain, deformation, and Okubo-Weiss
% parameter fields for the vector field with components U and V. The values
% dx and dy must be in meters, to scale final result to s^-1 (s^-2 for OWP).
% 
% Last Saved Time-stamp: <Fri 2008-12-05 14:26:26 Eastern Standard Time gramer>
%

  Nx = (size(U,2)-1) * dx;
  Ny = (size(U,1)-1) * dy;
  [X, Y] = meshgrid(0:dx:Nx, 0:dy:Ny);

  [CRLR,CAVR] = curl(X, Y, U, V);
  DIVR = divergence(X, Y, U, V);
  STRN = divergence(X, Y, U, -V);
  [STRS,CAVS] = curl(X, Y, -U, V);

  DFRM = sqrt( (STRN.^2) + (STRS.^2) );
  OWP = (DFRM.^2) - (CRLR.^2);

return;
