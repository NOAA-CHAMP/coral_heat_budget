function [T,Dt,B,Qv,Q] = thermal_exchange(t,s,h,q0,b,R,cfac,wfac)
% Function has now been superseded by HORIZONTAL_CONVECTION (vs.)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%function [T,Dt,B,Qv,Q] = thermal_exchange(t,s,h,q0,b,R,cfac,wfac)
%
% Thermally-induced exchange flow (thermal syphon, horizontal advection).
%
% From physical data (any combination of same-sized vectors and scalars):
%   t = sea temperature [oC]
%   s = salinity [psu]
%   h = site depth [m] or pressure [db]
%   q0 = net surface heat flux [W m^-2]
%   b = bottom slope (rise/run)
%   R = mixing(Ri) + entrainment rate (DEFAULT: 0.30)
%
% ... use formulae in literature to calculate thermal exchange parameters:
%   T = onset time for steady-state thermal exchange flow [s]
%   Dt = critical (maximum) depth for viscosity-dominated exchange [m]
%   B = net surface buoyance flux [m^2 s^-3]
%   Qv = volumetric discharge rate [m^2 s^-1]
%   Q = heat exchange [K/hr]
%
% SEE: Farrow and Patterson 1993, Sturman et al. 1999, Monismith et al. 2006,
%  Hughes and Griffiths 2008, Chubarenko 2010, Mao et al. 2009 and 2010
%
% Last Saved Time-stamp: <Mon 2011-04-18 14:02:10 Eastern Daylight Time gramer>

  if ( ~exist('R','var') || isempty(R) )
    R = 1 - 0.30;
  end;
  if ( ~exist('cfac','var') || isempty(cfac) )
    cfac = 1.20;
  end;
  if ( ~exist('wfac','var') || isempty(wfac) )
    wfac = 0.60;
  end;

  g = 9.79;
  alpha = sw_alpha(s,t,h);
  rho = sw_dens(s,t,h);
  cp = sw_cp(s,t,h);

  % Net surface buoyancy flux (ignoring runoff, for now...)
  B = (g .* alpha .* abs(q0)) ./ (rho .* cp);

  % Characteristic convective velocity scale
  uf = (B .* h) .^ (1/3);
  %DEBUG:
  disp({'M06 uf: ' nanmin(uf), nanmean(uf), nanmax(uf)});

  % Critical depth for viscous balance (assuming diurnal forcing)
  Dt = 0.3 .* uf .* (24*3600);


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Volumetric discharge rate Qv [m2/s]
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  % Sturman et al. 1999
  % l = 1000; %[m]
  % Qv = ((l .* B).^(1/3)) .* h;
  % Qv = 0.24 .* (B.^(1/3)) .* ( (l.*b./(1+b)) .^ (4/3) );
  % % Use h = l*b
  % Qv = 0.24 .* (B.^(1/3)) .* ( (h./(1+b)) .^ (4/3) );

  % Monismith et al. 2006
  % % Qv = 0.36 .* (1./(1+b)) .* uf .* h;
  % Qv = 0.3 .* uf .* h;
  % % Steady thermal balance, advective inertia
  % Qv = 0.15 .* (b.^(-1/3)) .* uf .* h;			% Poor
%   % Steady thermal balance, advective inertia
%   Qv = 1.0 .* (b.^(-1/3)) .* uf .* h;				% Poor

%   % Balanced thermal forcing, viscous/unsteady inertia
%   Qv = 1.00 .* sqrt((uf.^3) .* (24*3600) ./ h);		% Also not great

%   % Balanced thermal forcing, stress divergence
%   Qv = 1.00 .* uf;						% The worst


  % Unbalanced thermal forcing, viscous/unsteady inertia
  Qv = 1.00 .* b .* (uf.^3) .* ((24*3600)^2) ./ (h.^2);		% VERY GOOD for cooling

%   % Unbalanced thermal forcing, stress divergence
%   Qv = 1.00 .* sqrt(b .* (uf.^2) .* (24*3600) ./ h);		% OK for warming, terrible for cooling

%   % Unbalanced thermal forcing, advective inertia
%   %Qv = 1.00 .* sqrt((uf.^3) .* (24*3600) ./ h);		% NOTE TYPO in Monismith et al. 2006 Table 1!
%   Qv = 1.00 .* ((uf .* (24*3600) ./ h).^(3/2));		% Good for warming, OK for cooling


%   % Mix-and-match - but does *not* seem to perform better than "VERY GOOD" above
%   coolix = find(q0 < 0);
%   warmix = find(q0 > 0);
%   Qv = repmat(0,size(q0));
%   % Unbalanced thermal forcing, viscous/unsteady inertia
%   Qv(coolix) = 1.00 .* b .* (uf(coolix).^3) .* ((24*3600)^2) ./ (h(coolix).^2);		% VERY GOOD for cooling
%   % Unbalanced thermal forcing, advective inertia
%   Qv(warmix) = 1.0 .* ((uf(warmix) .* (24*3600) ./ h(warmix)).^(3/2));		% Good for warming, OK for cooling


%   % Use extremes as a proxy for onset lag times
%   bigix = find(abs(q0) > 50);
%   Qv = repmat(0,size(q0));
%   % Unbalanced thermal forcing, viscous/unsteady inertia
%   Qv(bigix) = 1.00 .* b .* (uf(bigix).^3) .* ((24*3600)^2) ./ (h(bigix).^2);		% VERY GOOD for cooling


%   % Mix-and-match - allow different flow (and moderation) rates for cooling vs. warming
%   coolix = find(q0 < 0);
%   warmix = find(q0 > 0);
%   Qv = repmat(0,size(q0));
%   % Unbalanced thermal forcing, viscous/unsteady inertia
%   Qv(coolix) = cfac .* b .* (uf(coolix).^3) .* ((24*3600)^2) ./ (h(coolix).^2);		% VERY GOOD for cooling
%   % Unbalanced thermal forcing, viscous/unsteady inertia
%   Qv(warmix) = wfac .* b .* (uf(warmix).^3) .* ((24*3600)^2) ./ (h(warmix).^2);		% VERY GOOD for cooling

  % Mix-and-match - allow different flow (and moderation) rates for cooling vs. warming
  coolix = find(q0 < 0);
  warmix = find(q0 > 0);
  % Unbalanced thermal forcing, viscous/unsteady inertia
  Qv(coolix) = cfac .* Qv(coolix);
  % Unbalanced thermal forcing, viscous/unsteady inertia
  Qv(warmix) = wfac .* Qv(warmix);

  %DEBUG:
  disp({'M06 Qv: ' nanmin(Qv), nanmean(Qv), nanmax(Qv)});


  % Onset time T for steady convection [s]

  % Chubarenko 2010
  l = 3600 .* Qv ./ h;
  D = h + (b.*l);
  dT = 3600 .* ( (q0./(cp.*rho.*(h+D))) - (q0./(cp.*rho.*h)) );
  % dRho = dT .* alpha,
  rho0 = sw_dens(s,(t+dT),h);
  dRho = (rho - rho0) ./ rho0;
  T = sqrt(l./(g.*abs(dRho).*abs(b)));
  %DEBUG:
  disp({'Chubarenko T: ' nanmin(T)/3600 nanmean(T)/3600 nanmax(T)/3600});


  % Lei and Patterson 2005
  nu = 1.05e-6;		% Molecular kinematic viscosity of seawater at 35psu, 20oC [m2/s]
  % nu = 1e-3;		% Eddy kinematic viscosity estimated over a reef [m2/s]
  k = 1.46e-7;		% Thermal diffusivity of seawater at 35psu, 20oC [m2/s]
  % Rac = 657.5;		% Critical Rayleigh number
  % Ra = g.*alpha.*(q0./(cp.*rho.*h)).*(h.^4)./(nu*(k^2));
  % T = sqrt(Rac/Ra).*(h^.2)./k;


  % Mao, Lei and Patterson 2009 (WARMING paper)
  Ra = g.*alpha.*(abs(q0)./(cp.*rho)).*(h.^4)./(nu*(k^2));
  %DEBUG:  disp({'MLP09 Ra: ' nanmin(Ra), nanmean(Ra), nanmax(Ra)});
  x = h ./ b;
  eta = 0.1;

  % T = (x.^(2/3)) .* (Ra.^(-1/3)) .* (l.^(4/3)) ./ k;
  T = (x.^(2/3)) .* exp(b.*eta.*x./3) .* (Ra.^(-1/3)) .* (h.^(4/3)) ./ k;
  %DEBUG:
  disp({'MLP09 T: ' nanmin(T)/3600 nanmean(T)/3600 nanmax(T)/3600});

  % Steady-state
  % NOTE (sloppy) TYPO in eqn. (22) of MLP10 WARMING paper, p. 180
  %u = ((x./h).^(1/3)) .* exp(-1./(3.*b.*x.*eta)) .* (Ra.^(1/3)) .* (k./h);
  % u = (x.^(1./3)).*exp(-b.*eta.*x./3).*(Ra.^(1/3)).*k.*(h.^(-4/3));
  % Unsteady
  T = (24*3600);
  u = Ra .* (T.^2) .* (k.^3) .* exp(-b.*eta.*x) ./ ( (h.^4) .* x );
  %DEBUG:
  disp({'MLP09 u: ' nanmin(u), nanmean(u), nanmax(u)});

  Dt_MLP09 = (x.^(1/3)) .* exp(b.*eta.*x./6) .* (h.^(2/3)) .* (Ra.^(-1/6));
  %DEBUG:
  disp({'MLP09 Dt: ' nanmin(Dt_MLP09), nanmean(Dt_MLP09), nanmax(Dt_MLP09)});

  %DEBUG:
  Dt_MLP09 = h;

  % Qv_MLP09 = (x.^(2/3)) .* exp(-b.*eta.*x./6) .* k .* (h.^(-2/3)) .* (Ra.^(1/6));
  Qv_MLP09 = u.*Dt_MLP09;
  %DEBUG:
  disp({'MLP09 Qv: ' nanmin(Qv_MLP09), nanmean(Qv_MLP09), nanmax(Qv_MLP09)});


%%%% ??? DEBUG
%   Qv(warmix) = wfac .* Qv_MLP09(warmix);
%%%% ??? DEBUG


  % Mao, Lei and Patterson 2010 (COOLING paper)
  %Ra = g.*alpha.*(abs(q0)./(cp.*rho.*h)).*(l.^4)./(nu*(k^2));
  Ra = g.*alpha.*(abs(q0)./(cp.*rho)).*(l.^4)./(nu*(k^2));
  x = h / b;
  % Time scale for reaching steady state of horizontal convection
  T = (x.^(2/3)) .* (Ra.^(-1/3)) .* (l.^(4/3)) ./ k;
  %DEBUG:
  disp({'MLP10 T: ' nanmin(T)/3600 nanmean(T)/3600 nanmax(T)/3600});

  % Transition distance from stable to unstable convective regime
  % lc = (Ra.^(-1/4)) .* (b.^(-3/2)) .* l;
  % %DEBUG:  [nanmin(lc), nanmean(lc), nanmax(lc)});

  % Steady state with an indistinct thermal boundary layer (SCALED h/k??)
  u = ( Ra .* (b.^4) .* (l.^-4) .* (x.^3) .* k ) .* (k./h);
  %DEBUG:
  disp({'MLP10 u: ' nanmin(u), nanmean(u), nanmax(u)});
  Qv_MLP10 = Ra .* (b.^5) .* ((x./l).^4) .* k .* (k./h);
  % Qv_MLP10 = u.*h;
  %DEBUG:
  disp({'MLP10 Qv: ' nanmin(Qv_MLP10), nanmean(Qv_MLP10), nanmax(Qv_MLP10)});

  % Steady state with a distinct thermal boundary layer (del_T < h)
  u = (Ra.^(1/3)) .* (l.^(-4/3)) .* (x.^(1/3)) .* k;
  %DEBUG:
  disp({'MLP10d u: ' nanmin(u), nanmean(u), nanmax(u)});
  Qv_MLP10 = (Ra.^(1/6)) .* ((x./l).^(2/3)) .* k;
  % This is almost certainly invalid!
  Qv_MLP10 = u.*h;
  %DEBUG:
  disp({'MLP10d Qv: ' nanmin(Qv_MLP10), nanmean(Qv_MLP10), nanmax(Qv_MLP10)});

  % Mao Lei Patterson u_f, Monismith Q_v
  Qv_MLP10nM = 1.00 .* b .* (u.^3) .* ((24*3600)^2) ./ (h.^2);		% VERY GOOD for cooling
  %DEBUG:  disp({'MLP10/M06 Qv: ' nanmin(Qv_MLP10nM), nanmean(Qv_MLP10nM), nanmax(Qv_MLP10nM)});


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Heat exchange rate Q [K/hr]
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  l = 3600 .* Qv ./ h;
  D = h + (b.*l);
  dT = 3600 .* ( (q0./(cp.*rho.*(h+D))) - (q0./(cp.*rho.*h)) );
  Q = R .* dT;

%%%% ??? DEBUG
  % Qv(T > (24*3600)) = 0;
  % Q(T > (24*3600)) = 0;
%%%% ??? DEBUG



return;
