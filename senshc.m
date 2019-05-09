1;
% Sensitivity analysis for Horizontal Convection

rhoCp = 4.09e6;
h = 3.55;
g = 9.8;
alph = 2.9e-4;

% for dt=[1 12]*3600;
% for doDt = [true false];
for dt=[12]*3600;
for doDt = [true];
% [bet,q0]=meshgrid(0.003:.001:0.20,0:10:1000);
[bet,q0]=meshgrid(0.003:.001:0.05,0:10:500);
B=(g.*alph.*abs(q0))./(rhoCp);
% Fischer et al (1979), Eq. 6.40 - per Monismith et al (2006)
uf=(B.*h).^(1/3);

% Steady thermal balance, advective inertia
Qv_SS=uf.*h./(bet.^(1/3));
% Balanced thermal forcing, viscous/unsteady inertia
Qv_SU=sqrt((uf.^3) .* dt ./ h);
% Unbalanced thermal forcing, advective inertia
Qv_US=(uf.*dt./h).^(3/2);
% Unbalanced thermal forcing, viscous/unsteady inertia
Qv_UU=bet.*(uf.^3).*(dt^2)./(h.^2);

for typ=1:4
switch (typ);
case 1; u=Qv_SS./h; ulbl='(\DeltaV_S_S)';
case 2; u=Qv_SU./h; ulbl='(\DeltaV_S_U)';
case 3; u=1.77e-4*Qv_US./h; ulbl='(\DeltaV_U_S)';
case 4; u=Qv_UU./h; ulbl='(\DeltaV_U_U)';
end;
dx=u*dt;
dTdtq0=dt.*q0./(rhoCp.*h);
dTdtx=dt.*q0./(rhoCp.*(h+(bet.*dx)));
dTdx=(dTdtq0-dTdtx)./dx;
dTdthc=-dt.*u.*dTdx;
dTdt=dTdtq0+dTdthc;

figure;
maxigraph;
ax(1) = gca;
if (doDt)
  [cf,hf]=contourf(q0,bet,dTdt.*3600./dt);
%  caxis([0 0.02]);
  titlename(['\partial_tT = (Q_0/\rhoC_ph) + u_H_C ^. \partial_xT_H_C vs. Q_0 and \beta ' ulbl]);
else
  [cf,hf]=contourf(q0,bet,dTdthc.*3600./dt);
  titlename(['\partial_tT_H_C = u_H_C ^. \partial_xT_H_C vs. Q_0 and \beta ' ulbl]);
end;
appendtitlename([' (after ' num2str(dt/3600) 'hrs)']);
cbh = colorbar('peer',ax(1)); ylabel(cbh,'\partial_tT [K/hr]');
ax(2) = axes('position',get(ax(1),'Position')); % STUPID colorbar behavior!
set(ax(2),'Color','none');
hold on;
linkaxes(ax);
[cu,hu]=contour(q0,bet,u,':','Color',[.5,.5,.5]); clabel(cu,hu,'Color',[.5,.5,.5]);
[cx,hx]=contour(q0,bet,dTdx,'k-.'); clabel(cx,hx,'Color','k');
xlabel('Q_0 [Wm^-^2]'); ylabel('\beta = \nablah');
lh = legend([hu,hx],'u_H_C [m/s]','\partial_xT_H_C [K/m]', 'Location','NorthEast');
set(lh,'Color','w');

end;
end;
end;

% figure;
% maxigraph;
% plot(q0,dTdtq0(:,1)*3600./dt);
% xlabel('Q_0 [Wm^-^2]'); ylabel('Q_0/\rhoC_ph [K/hr]');
