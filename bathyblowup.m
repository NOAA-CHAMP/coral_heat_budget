figspath = get_thesis_path('../figs');

for cst={'fwyf1','mlrf1','lonf1','smkf1','looe1','sanf1','dryf1'};
% for cst={'fwyf1','lonf1'};
  x=get_station_from_station_name(cst{:});
  x=get_ngdc_bathy_station(x);
  fmg;
  [ig,ig,ig,ig,cobj]=map_freef([x.lon-.4,x.lon+.4,x.lat-.4,x.lat+.4],'none');
  set(cobj,'LineWidth',1.75);
  [c,h]=contour(x.ngdc_92m_bathy.lon,x.ngdc_92m_bathy.lat,x.ngdc_92m_bathy.field,[-2,-10:-10:-100,-200:-100:-500],'k-','LineWidth',1);
  axis([x.lon-.15,x.lon+.15,x.lat-.14,x.lat+.14]);
  clabel(c,h);
  % contour(x.ngdc_92m_bathy.lon,x.ngdc_92m_bathy.lat,x.ngdc_92m_bathy.field,[-1:-1:-4],'k:');
  % contour(x.ngdc_92m_bathy.lon,x.ngdc_92m_bathy.lat,x.ngdc_92m_bathy.field,[-5:-10:-100],'k:','Color',[.5,.5,.5]);
  contour(x.ngdc_92m_bathy.lon,x.ngdc_92m_bathy.lat,x.ngdc_92m_bathy.field,[-5:-10:-100,-150:-100:-500],'k:','LineWidth',1,'LineColor',[.25,.25,.25]);
  plot(x.lon,x.lat,'kp','MarkerSize',14,'LineWidth',2,'MarkerFaceColor',[.5,.5,.5]);

  lm=axis;
  dlon=(abs(lm(1)-lm(2))/100);
  dlat=(abs(lm(3)-lm(4))/100);
  begpt=[lm(1),lm(3)]+[+dlon*10,+dlat*10];
  [endpt(1),endpt(2)]=transect_wgs84(begpt(1),begpt(2),1,90);
  midpt=mean([begpt;endpt]);
  lh(1)=line([begpt(1),endpt(1)],[begpt(2),endpt(2)]);
  lh(2)=line([begpt(1),begpt(1)],[begpt(2),begpt(2)+dlat]);
  lh(3)=line([endpt(1),endpt(1)],[endpt(2),endpt(2)+dlat]);
  set(lh,'LineWidth',2,'Color','k');
  text(midpt(1),midpt(2)-dlon,'1km','Color','k','FontSize',12,'FontWeight','bold','HorizontalAlignment','center');
  titlename(upper(cst{:}));
  print('-dtiff',fullfile(figspath,[cst{:},'-bathyblowup.tiff']));

  clear ig cobj c h lm dlon dlat begpt endpt midpt lh;
  x=[];
  clear x;
end;

clear cst figspath;
