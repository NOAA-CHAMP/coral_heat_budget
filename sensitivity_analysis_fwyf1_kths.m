doDiary=false;
doPrint=false;

diary off
more off;

figspath = get_thesis_path('../figs');

stnm='fwyf1';

RAPFX='erai';
KMPFX='avhrr_weekly';
ISPFX='ndbc';
TIDEPFX='tpxo_tide';
WAVEPFX='erai';
subs={};

% fname=[stnm,'-sensitivity-','control','.tif'],
% optimize_station_heat_budget(stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX);
% xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,fname)); end;

%for kth = {10,0,'ndbc_wind1_speed'};
for kth = {10};
  if ( ~ischar(kth{:}) )
    fname=[stnm,'-sensitivity-','ktheta',num2str(kth{:},'-%g'),'.tif'],
  else
    fname=[stnm,'-sensitivity-','ktheta','-',kth{:},'.tif'],
    stn = load_all_ndbc_data([],stnm);
    kth = { { stn.ndbc_wind1_speed,@(W)(min(10,((W./35).^2).*10)) } };
  end;
  optimize_station_heat_budget(stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,subs,struct('override_k_thetas',{kth}));
  xlabel(strrep(fname,'_','\_')); if (doPrint); print('-dtiff',fullfile(figspath,fname)); end;
  stn=[]; clear stn
end;
