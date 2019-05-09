1;
%% SCRIPT to draw regional-scale land-bathymetry map of southern Florida
%
% Last Saved Time-stamp: <Fri 2012-03-02 15:11:23  lew.gramer>


if ( ~exist('fwyf1','var') )
  lkwf1=get_station_from_station_name('lkwf1');
  fwyf1=get_station_from_station_name('fwyf1');
  mlrf1=get_station_from_station_name('mlrf1');
  lonf1=get_station_from_station_name('lonf1');

  % SE Straits of Florida - for isobaths south of MLRF1
  sesf1.station_name='sesf1'; sesf1.lon=-80.5; sesf1.lat=24.3;

  % SW Florida Shelf - for isobaths north & west of NFBF1
  swfs1.station_name='swfs1'; swfs1.lon=-81.6; swfs1.lat=25.3;

  % Central W Florida Shelf - for isobaths north of SWFS1
  cwfs1.station_name='cwfs1'; cwfs1.lon=-81.6; cwfs1.lat=26.0;

  % Offshore of the Marquesas - for isobaths btw SANF1 and DRYF1
  marf1.station_name='marf1'; marf1.lon=-82.40; marf1.lat=24.57;

  nfbf1=get_station_from_station_name('nfbf1');
  smkf1=get_station_from_station_name('smkf1');
  looe1=get_station_from_station_name('looe1');
  sanf1=get_station_from_station_name('sanf1');
  dryf1=get_station_from_station_name('dryf1');
end;

% % isodepths = [-2,-6,-10,-30,-80,-200,-300,-500,-700];
% isodepths = [-2,-10,-30,-80,-200,-300,-700];
isodepths = [-2,-10,-30,-80,-150,-300,-700];
isocolor = [.4,.4,.4];


fmg;

%% Draw land outline

% % For ICRS paper
% %map_freef([mlrf1.lon-0.97,mlrf1.lon+0.97,mlrf1.lat-1.00,mlrf1.lat+1.00],isodepths);
% [ig,ig,cs,ch,cobj,clhs]=map_freef([stn.lon-0.66,stn.lon+0.66,stn.lat-0.60,stn.lat+0.60],isodepths);

% For Gramer & Mariano
%[ig,ig,cs,ch,cobj,clhs]=map_freef([-82.1,-79.9,24.2,25.8],'none');
[ig,ig,cs,ch,cobj,clhs]=map_freef([-82.9,-79.9,24.2,25.8],'none');

set(cobj,'LineWidth',1.5);


%% Draw NGDC isobaths for 40x40km square around each site

if ( ~isfield(fwyf1,'ngdc_92m_bathy') ); fwyf1 = get_ngdc_bathy_station(fwyf1); end;
[cs,ch] = contour(fwyf1.ngdc_92m_bathy.lon, fwyf1.ngdc_92m_bathy.lat, fwyf1.ngdc_92m_bathy.field, ...
                  isodepths, 'LineColor', isocolor); clhs = clabel(cs,ch,'Color',isocolor,'FontSize',6);

if ( ~isfield(mlrf1,'ngdc_92m_bathy') ); mlrf1 = get_ngdc_bathy_station(mlrf1); end;
[cs,ch] = contour(mlrf1.ngdc_92m_bathy.lon, mlrf1.ngdc_92m_bathy.lat, mlrf1.ngdc_92m_bathy.field, ...
                  isodepths, 'LineColor', isocolor); %clhs = clabel(cs,ch,'Color',isocolor,'FontSize',6);

if ( ~isfield(lonf1,'ngdc_92m_bathy') ); lonf1 = get_ngdc_bathy_station(lonf1); end;
[cs,ch] = contour(lonf1.ngdc_92m_bathy.lon, lonf1.ngdc_92m_bathy.lat, lonf1.ngdc_92m_bathy.field, ...
                  isodepths, 'LineColor', isocolor); clhs = clabel(cs,ch,'Color',isocolor,'FontSize',6);

if ( ~isfield(nfbf1,'ngdc_92m_bathy') ); nfbf1 = get_ngdc_bathy_station(nfbf1); end;
[cs,ch] = contour(nfbf1.ngdc_92m_bathy.lon, nfbf1.ngdc_92m_bathy.lat, nfbf1.ngdc_92m_bathy.field, ...
                  isodepths, 'LineColor', isocolor); clhs = clabel(cs,ch,'Color',isocolor,'FontSize',6);

if ( ~isfield(sesf1,'ngdc_92m_bathy') ); sesf1 = get_ngdc_bathy_station(sesf1); end;
[cs,ch] = contour(sesf1.ngdc_92m_bathy.lon, sesf1.ngdc_92m_bathy.lat, sesf1.ngdc_92m_bathy.field, ...
                  isodepths, 'LineColor', isocolor); clhs = clabel(cs,ch,'Color',isocolor,'FontSize',6);

if ( ~isfield(swfs1,'ngdc_92m_bathy') ); swfs1 = get_ngdc_bathy_station(swfs1); end;
[cs,ch] = contour(swfs1.ngdc_92m_bathy.lon, swfs1.ngdc_92m_bathy.lat, swfs1.ngdc_92m_bathy.field, ...
                  isodepths, 'LineColor', isocolor); clhs = clabel(cs,ch,'Color',isocolor,'FontSize',6);

if ( ~isfield(cwfs1,'ngdc_92m_bathy') ); cwfs1 = get_ngdc_bathy_station(cwfs1); end;
[cs,ch] = contour(cwfs1.ngdc_92m_bathy.lon, cwfs1.ngdc_92m_bathy.lat, cwfs1.ngdc_92m_bathy.field, ...
                  isodepths, 'LineColor', isocolor); %clhs = clabel(cs,ch,'Color',isocolor,'FontSize',6);

if ( ~isfield(smkf1,'ngdc_92m_bathy') ); smkf1 = get_ngdc_bathy_station(smkf1); end;
[cs,ch] = contour(smkf1.ngdc_92m_bathy.lon, smkf1.ngdc_92m_bathy.lat, smkf1.ngdc_92m_bathy.field, ...
                  isodepths, 'LineColor', isocolor); %clhs = clabel(cs,ch,'Color',isocolor,'FontSize',6);

if ( ~isfield(looe1,'ngdc_92m_bathy') ); looe1 = get_ngdc_bathy_station(looe1); end;
[cs,ch] = contour(looe1.ngdc_92m_bathy.lon, looe1.ngdc_92m_bathy.lat, looe1.ngdc_92m_bathy.field, ...
                  isodepths, 'LineColor', isocolor); %clhs = clabel(cs,ch,'Color',isocolor,'FontSize',6);

if ( ~isfield(sanf1,'ngdc_92m_bathy') ); sanf1 = get_ngdc_bathy_station(sanf1); end;
[cs,ch] = contour(sanf1.ngdc_92m_bathy.lon, sanf1.ngdc_92m_bathy.lat, sanf1.ngdc_92m_bathy.field, ...
                  isodepths, 'LineColor', isocolor); clhs = clabel(cs,ch,'Color',isocolor,'FontSize',6);

if ( ~isfield(marf1,'ngdc_92m_bathy') ); marf1 = get_ngdc_bathy_station(marf1); end;
[cs,ch] = contour(marf1.ngdc_92m_bathy.lon, marf1.ngdc_92m_bathy.lat, marf1.ngdc_92m_bathy.field, ...
                  isodepths, 'LineColor', isocolor); %clhs = clabel(cs,ch,'Color',isocolor,'FontSize',6);

if ( ~isfield(dryf1,'ngdc_92m_bathy') ); dryf1 = get_ngdc_bathy_station(dryf1); end;
[cs,ch] = contour(dryf1.ngdc_92m_bathy.lon, dryf1.ngdc_92m_bathy.lat, dryf1.ngdc_92m_bathy.field, ...
                  isodepths, 'LineColor', isocolor); %clhs = clabel(cs,ch,'Color',isocolor,'FontSize',6);


daspect([1/cosd(25),1,1]);

%text(-80.8,25.5,'Florida', 'FontSize',24, 'FontWeight','bold');
text(-80.9,25.6,'Florida', 'FontSize',24, 'FontWeight','normal');

plot(fwyf1.lon,fwyf1.lat,'kp', 'MarkerSize',10, 'MarkerFaceColor','k');
plot(mlrf1.lon,mlrf1.lat,'kp', 'MarkerSize',10, 'MarkerFaceColor','k');
plot(lonf1.lon,lonf1.lat,'kp', 'MarkerSize',10, 'MarkerFaceColor','k');
plot(smkf1.lon,smkf1.lat,'kp', 'MarkerSize',10, 'MarkerFaceColor','k');
plot(looe1.lon,looe1.lat,'kp', 'MarkerSize',10, 'MarkerFaceColor','k');
plot(sanf1.lon,sanf1.lat,'kp', 'MarkerSize',10, 'MarkerFaceColor','k');
plot(dryf1.lon,dryf1.lat,'kp', 'MarkerSize',10, 'MarkerFaceColor','k');

text(fwyf1.lon,fwyf1.lat,'. ^\leftarrow FWYF1','Rotation',+000, 'FontSize',8);
text(mlrf1.lon,mlrf1.lat,'. ^\leftarrow MLRF1','Rotation',-045, 'FontSize',8);
text(lonf1.lon,lonf1.lat,'. ^\leftarrow LONF1','Rotation',+090, 'FontSize',8);
text(smkf1.lon,smkf1.lat,'. ^\leftarrow SMKF1','Rotation',-060, 'FontSize',8);
text(looe1.lon,looe1.lat,'. ^\leftarrow LOOE1','Rotation',-068, 'FontSize',8);
text(sanf1.lon,sanf1.lat,'. ^\leftarrow SANF1','Rotation',-079, 'FontSize',8);
text(dryf1.lon,dryf1.lat,'. ^\leftarrow DRYF1','Rotation',-055, 'FontSize',8);

if (1)
  print('-dtiff',fullfile(get_thesis_path('../figs'),'Gramer & Mariano - Fig1.tiff'));
else
  disp('DID *NOT* PRINT THIS FIGURE YET...');
end;
