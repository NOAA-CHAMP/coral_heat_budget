1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SCRIPT to analyze AOML South Florida Program Thermosalinograph data for
% cross-shore gradients in near-surface sea temperature, in particular
% relative to University of South Florida 1km AVHRR weekly fields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% SCRIPT used to calculate values for Appendix II of Lew Gramer's Dissertation
%
% Last Saved Time-stamp: <Mon 2013-07-01 20:09:17 Eastern Daylight Time gramer>


tic,
more off;

figspath = get_thesis_path('../figs');
DISSpath = get_thesis_path('../DISS');

%stnm = 'fwyf1';
stnm = 'mlrf1';
%stnm = 'smkf1';
%stnm = 'looe1';
%stnm = 'sanf1';
%stnm = 'dryf1';

% % for cstnm = {
% %     'fwyf1' ...
% %     'bnppa' ...
% %     'cryf1' ...
% %     'mlrf1' ...
% %     'tnrf1' ...
% %     'smkf1' ...
% %     'looe1' ...
% %     'amsf1' ...
% %     'sanf1' ...
% %     'plsf1' ...
% %             };
% for cstnm = {
%     'fwyf1' ...
%     'mlrf1' ...
%     'smkf1' ...
%     'looe1' ...
%     'sanf1' ...
%             };
% stnm = cstnm{:};


%%%%%%%%%%
%% Load data for analysis

if ( ~exist('stn','var') || ~isfield(stn,'station_name') || ~strcmpi(stn.station_name,stnm) )
  stn = []; clear stn;
  stn = get_station_from_station_name(stnm);
  stn = station_optimal_isobath_orientation(stn);
  if ( strcmpi(stnm,'sanf1') )
    stn = get_avhrr_weekly_field(stn,true,'nearest',3);
  elseif ( strcmpi(stnm,'dryf1') )
    stn = get_avhrr_weekly_field(stn,true,{@nanmean,2,2},3);
  else
    stn = get_avhrr_weekly_field(stn,true,'linear',5);
  end;
  stn = station_reorient_vectors(stn,stn.isobath_orientation,...
                                 'hourly_avhrr_weekly_sst_x',...
                                 'hourly_avhrr_weekly_sst_y');
end;

if ( ~exist('tsg','var') )
    tsg = antsg;
end;


%%%%%%%%%%
%% Selection criteria for shelf-break, onshore, and offshore points

% Find all TSG data points at the "shelf-break" (20m isobath)
shelf_break_depth = -20;
brkix = find(roundn(tsg.z(:),0) == shelf_break_depth);

% Find all TSG "shelf-break" data points near station STNM (<25km)
%minstndx = 25;
minstndx = 30;
[rng,az]=distance_wgs84(stn.lat,stn.lon,tsg.lat(brkix),tsg.lon(brkix));
brkix = brkix(abs(rng) <= minstndx);


% Find all TSG points approx. 1km inshore or offshore of the shelf-break

% 144 secs steaming at top speed directly inshore should take you <=1km
%mindt = (144/3600/24);
mindt = (20/3600/24);
% After 6 hours, we may argue cross-shore sea temperature gradients decorrelate
%maxdt = (6/24);
maxdt = 1;

% We want points ~1km in- and ~1km offshore of the shelf-break
% mindx = 0.8; maxdx = 1.2;
mindx = 0.5; maxdx = 1.5;

% "Inshore" means less than ~16m depth, "offshore" < ~26m depth.
% % minz = -10; maxz = -40;
% minz = shelf_break_depth + 6;
% maxz = shelf_break_depth - 6;
minz = shelf_break_depth + 4;
maxz = shelf_break_depth - 4;

% Also use relative azimuthal heading from shelf-break to check "in-" and
% "offshore": heading should be perpendicular to isobath +/- 15 degrees
%min_rel_hdg = cosd(30);
min_rel_hdg = cosd(50);



%%%%%%%%%%
%% Analyze TSG gradients near station relative to AVHRR weekly SST data

missed = 0;
brksst.date = tsg.date(brkix);
brksst.data = repmat(nan,[numel(brkix) 1]);
onsst.date = tsg.date(brkix);
onsst.data = repmat(nan,[numel(brkix) 1]);
onsst.disp = repmat(nan,[numel(brkix) 1]);
offsst.date = tsg.date(brkix);
offsst.data = repmat(nan,[numel(brkix) 1]);
offsst.disp = repmat(nan,[numel(brkix) 1]);

all_onsst.date = []; all_onsst.data = []; all_onsst.lon = []; all_onsst.lat = []; 
all_offsst.date = []; all_offsst.data = []; all_offsst.lon = []; all_offsst.lat = []; 

dxssst.date = tsg.date(brkix);
dxssst.data = repmat(nan,[numel(brkix) 1]);
dxssst.disp = repmat(nan,[numel(brkix) 1]);

brklons = []; brklats = [];
onlons = []; onlats = []; onzs = [];
offlons = []; offlats = []; offzs = [];

kept_brkix = brkix;
for ix = 1:length(brkix)
    dt = abs(tsg.date-tsg.date(brkix(ix)));
    %nrix = find((10/3600/24)<=dt & dt<=1);
    nronix = find(mindt<=dt & dt<=maxdt & tsg.z<=-1 & tsg.z>=minz);
    nroffix = find(mindt<=dt & dt<=maxdt & tsg.z<=-1 & tsg.z<=maxz);
    dt=[]; clear dt;

    [rng,az]=distance_wgs84(tsg.lat(brkix(ix)),tsg.lon(brkix(ix)),tsg.lat(nronix),tsg.lon(nronix));
    relhdg = cosd( (stn.isobath_orientation-90) - az );
    nronix = nronix(mindx<=abs(rng) & abs(rng)<=maxdx & relhdg>=min_rel_hdg);

    [rng,az]=distance_wgs84(tsg.lat(brkix(ix)),tsg.lon(brkix(ix)),tsg.lat(nroffix),tsg.lon(nroffix));
    relhdg = cosd( (stn.isobath_orientation+90) - az );
    nroffix = nroffix(mindx<=abs(rng) & abs(rng)<=maxdx & relhdg>=min_rel_hdg);

    if ( isempty(nronix) || isempty(nroffix) )
        missed = missed + 1;
        kept_brkix(kept_brkix==brkix(ix)) = [];
    else
        brksst.data(ix) = tsg.sst(brkix(ix));
        onsst.data(ix) = nanmean(tsg.sst(nronix));
        onsst.disp(ix) = nanstd(tsg.sst(nronix));
        offsst.data(ix) = nanmean(tsg.sst(nroffix));
        offsst.disp(ix) = nanstd(tsg.sst(nroffix));

        all_onsst.date(end+1:end+length(nronix)) = tsg.date(brkix(ix));
        all_onsst.data(end+1:end+length(nronix)) = tsg.sst(nronix);
        all_onsst.lon(end+1:end+length(nronix)) = tsg.lon(nronix);
        all_onsst.lat(end+1:end+length(nronix)) = tsg.lat(nronix);
        all_offsst.date(end+1:end+length(nroffix)) = tsg.date(brkix(ix));
        all_offsst.data(end+1:end+length(nroffix)) = tsg.sst(nroffix);
        all_offsst.lon(end+1:end+length(nroffix)) = tsg.lon(nroffix);
        all_offsst.lat(end+1:end+length(nroffix)) = tsg.lat(nroffix);

        % Cross-shore SST gradient in [K/m]
        dxssst.data(ix) = (offsst.data(ix) - onsst.data(ix)) / 2e3;
        %dxssst.data(ix) = (offsst.data(ix) - brksst.data(ix)) / 1e3;
        dxssst.disp(ix) = sqrt((offsst.disp(ix).^2) - (onsst.data(ix).^2)) / 2e3;

        % Diagnostic stuff
        brklons = [brklons(:) ; tsg.lon(brkix(ix))];
        brklats = [brklats(:) ; tsg.lat(brkix(ix))];
        onlons = [onlons(:) ; tsg.lon(nronix)];
        onlats = [onlats(:) ; tsg.lat(nronix)];
        onzs = [onzs(:) ; tsg.z(nronix)];
        offlons = [offlons(:) ; tsg.lon(nroffix)];
        offlats = [offlats(:) ; tsg.lat(nroffix)];
        offzs = [offzs(:) ; tsg.z(nroffix)];
    end;
end;

missed,


%%%%%%%%%%
%% Plot results

if ( ~isempty(onlons) )
    if (0)
      fmg; plot_ts(subset_ts(dxssst,@ts_boreal_cool),'b.',subset_ts(dxssst,@ts_boreal_warm),'r.');
      ylim([-1e-3,+1e-3]);
      titlename([upper(stn.station_name),' TSG sea-temperature gradient time-series']);
      fmg; grpplot_ts(dxssst,@get_week,@nanmean,0,'*');
      ylim([-1e-3,+1e-3]);
      titlename([upper(stn.station_name),' TSG sea-temperature gradient climatology (N=',num2str(length(dxssst.data)),' points)']);
    end;

    if (1-1)
      fmg; boxplot_ts(dxssst);
      ylim([-1e-3,+1e-3]);
      titlename([upper(stn.station_name),' TSG sea-temperature gradient climatology']);
    end;

    if (0)
      fmg; plot_ts(brksst,'k.',all_onsst,'r.',all_offsst,'b.'); datetick3('x',27,'keeplimits');
      titlename([upper(stn.station_name),' TSG sea temperatures']); legend('Shelf-break','Inshore','Offshore');
      fmg; plot_ts(brksst,'k.',onsst,'r.',offsst,'b.'); datetick3('x',27,'keeplimits');
      titlename([upper(stn.station_name),' TSG mean sea temperature']); legend('Shelf-break','Inshore','Offshore');
    end;

    if (1-1+0)
      stn = get_ngdc_bathy_station(stn);
      stn = plot_ngdc_bathy_station(stn,-[2:2:10,15:5:40],[],true,@contour);
      %plot(brklons,brklats,'k.'); plot(onlons,onlats,'r.'); plot(offlons,offlats,'b.');
      plot(brklons,brklats,'k.');
      plot3(all_onsst.lon,all_onsst.lat,get_month(all_onsst.date),'r.');
      plot3(all_offsst.lon,all_offsst.lat,get_month(all_offsst.date),'b.');
      plot(stn.lon,stn.lat,'kp');
      set_pcolor_cursor;
      titlename([upper(stn.station_name),' TSG data points ',num2str(get_year(tsg.date([1,end]))')]);
    end;

    if (1-1)
      % scatter_fit_ts_seasons(stn.hourly_avhrr_weekly_sst_xshore,dxssst,[],[],'AVHRR Weekly','\partial_x_sTSG',[],[],true);
      % subplots_set('xlim',[-1e-3,+1e-3],'ylim',[-1e-3,+1e-3]); appendtitlename([' ',upper(stn.station_name)]);
      scatter_fit_ts(stn.hourly_avhrr_weekly_sst_xshore,dxssst,[],@ts_boreal_warm,'AVHRR Weekly','\partial_x_sTSG',[],[],true);
      axis([-1e-3,+1e-3,-1e-3,+1e-3]); appendtitlename([' ',upper(stn.station_name)]);
      print('-dtiff',fullfile(figspath,[stnm,'-avhrr_weekly_sst_xshore-scatter-SFP_TSG-warm.tif']));

      scatter_fit_ts(stn.hourly_avhrr_weekly_sst_xshore,dxssst,[],@ts_boreal_cool,'AVHRR Weekly','\partial_x_sTSG',[],[],true);
      axis([-1e-3,+1e-3,-1e-3,+1e-3]); appendtitlename([' ',upper(stn.station_name)]);
      print('-dtiff',fullfile(figspath,[stnm,'-avhrr_weekly_sst_xshore-scatter-SFP_TSG-cool.tif']));
    end;

    if (1)
      scatter_fit_ts(stn.hourly_avhrr_weekly_sst_xshore,dxssst,[],[],'AVHRR Weekly','\partial_x_sTSG',[],[],true);
      print('-dtiff',fullfile(DISSpath,[stnm,'-avhrr_weekly_sst_xshore-scatter-SFP_TSG.tif']));

      scatter_fit_ts_seasons(stn.hourly_avhrr_weekly_sst_xshore,dxssst,[],[],'AVHRR Weekly','\partial_x_sTSG',[],[],true);
      subplots_set('xlim',[-1e-3,+1e-3],'ylim',[-1e-3,+1e-3]); appendtitlename([' ',upper(stn.station_name)]);
      print('-dtiff',fullfile(DISSpath,[stnm,'-avhrr_weekly_sst_xshore-scatter-SFP_TSG-seasons.tif']));
    end;

end;

pause(1);

% end; %for cstnm

% clear DISSpath all_o* az brk* cfn dx* figspath fn grid_interp_method ix kept_brkix
% clear maxd* maxz min_rel* mind* minz missed npts nro* offl* offsst minstndx offz*
% clear onlats onlons onsst onzs relhdg rng shelf_break_depth

toc,
more on;
