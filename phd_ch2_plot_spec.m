1;

tic,
more off;

printFmt='tiff';
%printFmt=[];

doNorm = true;

%xlm=[];
xlm=[1e-1,3e3];

def_ylm=[1e-3,1e6];
def_ylm_ts=[1e-5,1e6];
def_ylm_td=[1e-5,1e4];
def_ylm_cu=[1e-5,1e2];
def_ylm_qa=[1e-9,1e2];

%def_lineWidth=1.5;
def_lineWidth=2.5;

%for cnm={'fwyf1','mlrf1','lonf1','smkf1','sanf1','dryf1'};
for cnm={'lonf1','smkf1'};
%for cnm={'smkf1'};
    nm=cnm{:};
    stn=get_station_from_station_name(nm);
    stn=load_all_ndbc_data(stn);
    if (isfield(stn,'ndbc_dew_t'));
        stn=station_dewp_to_relhumid(stn,'ndbc_air_t','ndbc_dew_t','ndbc_relhumid');
        stn=station_relhumid_to_spechumid(stn,'ndbc_air_t','ndbc_relhumid','ndbc_spechumid');
    end;

    stn=station_optimal_isobath_orientation(stn);
    stn=verify_variable(stn,{'ndbc_wind1_u','ndbc_wind1_v'});
    stn=station_reorient_vectors(stn,'isobath_orientation','ndbc_wind1_u','ndbc_wind1_v');

    %for cfld={'ndbc_sea_t','ndbc_air_t','ndbc_barom','ndbc_wind1_speed','ndbc_wind1_xshore','ndbc_wind1_lshore','ndbc_tide','ndbc_spechumid'};
    for cfld={'ndbc_tide'};
        fld=cfld{:};
        if (isfield(stn,fld));
            disp([nm,'.',fld]);
            ylm=def_ylm;
            lineWidth=def_lineWidth;
            switch (lower(fld)),
             case 'ndbc_spechumid',
              ylm=def_ylm_qa;
             case 'ndbc_tide',
              ylm=def_ylm_td;
             case 'ndbc_sea_t',
              ylm=def_ylm_ts;
              % lineWidth=2.5;
            end;
            [Pxx,Wd,fh,lhs]=plot_spec(stn,fld,[],[],xlm,ylm,printFmt,doNorm);
            set(lhs,'LineWidth',lineWidth);
        end;
    end;
    %%%%???DEBUG
    if ( ~isempty(printFmt) ); close all; end;
    stn=[]; clear stn Pxx Wd fh lhs;
end;

if (0)

ylm=def_ylm;
lineWidth=def_lineWidth;

stn = get_station_from_station_name('looe1');
stn = get_looe1_microcat(stn);
stn = get_looe1_adcp(stn);
% Contiguous (gap<=1d) ocean currents indices
contigix=1959:30096;
for cfld={'microcat_seatemp','adcp_seatemp','adcp_speed','adcp_sfc_speed','adcp_x','adcp_l'};
    fld=cfld{:};
    switch (lower(fld)),
     case 'microcat_seatemp',
      ylm=def_ylm_ts;
     case 'adcp_seatemp',
      ylm=def_ylm_ts;
     otherwise,     % CURRENTS
      ylm=def_ylm_cu;
      stn.(['contig_',fld]) = subset_ts(stn.(fld),contigix);
      fld=['contig_',fld];
    end;
    [Pxx,Wd,fh,lhs]=plot_spec(stn,fld,[],[],xlm,ylm,printFmt,doNorm);
    set(lhs,'LineWidth',lineWidth);
end;
%%%%???DEBUG
if ( ~isempty(printFmt) ); close all; end;
stn=[]; clear stn Pxx Wd fh lhs;

end;

clear cnm nm
more on;
toc,
