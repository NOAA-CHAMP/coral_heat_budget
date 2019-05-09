1;
th=0:0.01:6;
f=sw_f(25);
Ts=[24.9 24.5:-0.5:15];
for nix=1:numel(Ts); Ns(nix)=sqrt(sw_bfrq([35,35],[25,Ts(nix)],[5,30])); end;
fmg; cm=colormap;
for nix=1:numel(Ns);
    N=Ns(nix);
    plot(th,(2*pi/3600/24)./sqrt(((N.*sind(th)).^2)+((f.*cosd(th)).^2)),'Color',cm(nix*2,:));
end;
xlabel('\theta'); ylabel('2\pi/\omega_c [days]'); legend(num2str(Ns','N=%g'));

