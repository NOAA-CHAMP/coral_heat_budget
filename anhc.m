1;

  s = 35;
  t = 25;
  h = 10.96;
  dt = 3600;

  R = (1-0.08);

  rho = sw_dens(s,t,h);					%[kg/m^3]
  Cp = sw_cp(s,t,h);					%[J/kg*K]
  rhoCp = rho.*Cp;						%[J/K*m^3]
  g = 9.79;								%[m/s^2]
  alph = sw_alpha(s,t,h);				%[1/K]

  q0s = linspace(-500,0,20);
  %q0s = linspace(-200,0,20);
  %bets = [0.002,0.02,0.20];
  bets = 2.*logspace(-3,-1,20);

  clear dTdthc_SS dTdthc_US dTdthc_SU dTdthc_UU
  for qix=1:numel(q0s);
      q0 = q0s(qix);
      for bix=1:numel(bets);
          bet = bets(bix);
          res.B0 = (g.*alph.*abs(q0))./(rhoCp);
          res.uf = (res.B0.*h).^(1/3);
          res.Qv_SS = res.uf.*h./(bet.^(1/3));
          res.Qv_US = sqrt((res.uf.^3) .* (24.*dt) .* h);
          res.Qv_SU = (bet.^(2/3)).*res.uf.*((res.uf.*(24.*dt)./h).^(3/2));
          res.Qv_UU = bet.*(res.uf.^3).*((24.*dt)^2)./h;

          res.u_SS  = (5.0.*res.Qv_SS./h) - 0.05;   res.u_SS(res.u_SS<0) = 0;
          res.u_US  = (3.0.*res.Qv_US./h) - 0.026;  res.u_US(res.u_US<0) = 0;
          res.u_SU  = (2.7.*res.Qv_SU./h) - 0.0224; res.u_SU(res.u_SU<0) = 0;
          res.u_UU  = (0.1.*res.Qv_UU./h);          res.u_UU(res.u_UU<0) = 0;

          res.dTdtq0 = dt.*q0./(rhoCp.*h);

          res.dx_SS = res.u_SS .* dt;
          res.dx_US = res.u_US .* dt;
          res.dx_SU = res.u_SU .* dt;
          res.dx_UU = res.u_UU .* dt;

          res.dTdtx_SS = dt.*q0./(rhoCp.*(h+(bet.*res.dx_SS)));
          res.dTdtx_US = dt.*q0./(rhoCp.*(h+(bet.*res.dx_US)));
          res.dTdtx_SU = dt.*q0./(rhoCp.*(h+(bet.*res.dx_SU)));
          res.dTdtx_UU = dt.*q0./(rhoCp.*(h+(bet.*res.dx_UU)));

          res.dTdx_SS=(res.dTdtq0-res.dTdtx_SS)./res.dx_SS;
          res.dTdx_SS(~isfinite(res.dTdx_SS)) = 0;
          res.dTdx_US=(res.dTdtq0-res.dTdtx_US)./res.dx_US;
          res.dTdx_US(~isfinite(res.dTdx_US)) = 0;
          res.dTdx_SU=(res.dTdtq0-res.dTdtx_SU)./res.dx_SU;
          res.dTdx_SU(~isfinite(res.dTdx_SU)) = 0;
          res.dTdx_UU=(res.dTdtq0-res.dTdtx_UU)./res.dx_UU;
          res.dTdx_UU(~isfinite(res.dTdx_UU)) = 0;

          dTdthc_SS(bix,qix) = -R.*24.*dt.*res.u_SS.*res.dTdx_SS;
          dTdthc_US(bix,qix) = -R.*24.*dt.*res.u_US.*res.dTdx_US;
          dTdthc_SU(bix,qix) = -R.*24.*dt.*res.u_SU.*res.dTdx_SU;
          dTdthc_UU(bix,qix) = -R.*24.*dt.*res.u_UU.*res.dTdx_UU;
      end;
  end;

  fh=fmg;
  surf(q0s,bets,dTdthc_SS);
  surf(q0s,bets,dTdthc_US);
  surf(q0s,bets,dTdthc_SU);
  surf(q0s,bets,dTdthc_UU);
  set(gca,'yscale','log');
  xlim([-500,0]);
  ylim([0.002,0.2]);
  zlim([0,1]);
  %view(285,30);
  view(618,-324);
  shading interp;

  text(q0s(1),bets(end)*1.0,dTdthc_SS(end,1),'SS');
  text(q0s(1),bets(end)*1.2,dTdthc_US(end,1),'US');
  text(q0s(1),bets(end)*1.4,dTdthc_SU(end,1),'SU');
  text(q0s(1),bets(end)*1.6,dTdthc_UU(end,1),'UU');

  xlabel('Surface cooling [Wm^-^2]');
  ylabel('Bottom slope \beta');
  zlabel('Horizontal convection [Kd^-^1]');


  set(gca,'zlim',[0,1]);

  titlename('Sensitivity analysis: Horizontal convection vs. \beta, Q_0');

  print('-dtiff',fullfile(get_thesis_path('../figs'),'phd-hc-vs-beta-vs-q0.tif'));
  saveas(fh,fullfile(get_thesis_path('../figs'),'phd-hc-vs-beta-vs-q0.fig'));

  clear s t h dt R rho Cp rhoCp g alph qix bix q0 bet res ig zx ix jx fh ans
  %clear q0s bets
