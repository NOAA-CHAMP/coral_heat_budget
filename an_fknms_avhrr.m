1;

if ( ~exist('stns','var') )
    stns = get_fknms_thermistors;
end;

for cst=grepstruct(stns,'FKNMS')';
    st=cst{:}; 
    if (isfield(stns.(st),'fknms_seatemp') && numel(stns.(st).fknms_seatemp.data)>=(12*365));
        if ( ~isfield(stns.(st),'raw_avhrr_weekly_sst') )
            stns.(st) = get_avhrr_weekly_field(stns.(st));
            stns.(st) = rmfield(stns.(st),grepstruct(stns.(st),'_field'));
        end;
        disp(st);
        stns.(st)=verify_variable(stns.(st),{'fknms_seatemp_1_d_avg','fknms_seatemp_7_d_avg'},true);
        stns.(st).fknms_seatemp_1_d_avg.date = stns.(st).fknms_seatemp_1_d_avg.date - 0.5;
        stns.(st).fknms_seatemp_7_d_avg.date = stns.(st).fknms_seatemp_7_d_avg.date - 3.5;
        try,
            fh=fmg; set(gca,'FontSize',6);
            %scatter_fit_ts(stns.(st).raw_avhrr_weekly_sst,...
            scatter_fit_ts_seasons(stns.(st).raw_avhrr_weekly_sst,...
                stns.(st).fknms_seatemp_1_d_avg,[],[],'AVHRR',...
                ['\mu_7_d',strrep(upper(st),'_','\_'),'-3.5d'],fh,[],true);
            axis([14,34,14,34]);
            subplots_set(fh,'FontSize',6,'XLim',[14,34],'YLim',[14,34]);
            %print('-dtiff',['../figs/',lower(st),'-avhrr_weekly_sst-scatter-seasons-fknms_seatemp_1_d_avg.tif']);
            %%print('-dtiff',['../figs/',lower(st),'-avhrr_weekly_sst-scatter-fknms_seatemp_1_d_avg.tif']);
            print('-dtiff',['../figs/',lower(st),'-avhrr_weekly_sst-scatter-seasons-fknms_seatemp_7_d_avg.tif']);
            %print('-dtiff',['../figs/',lower(st),'-avhrr_weekly_sst-scatter-fknms_seatemp_7_d_avg.tif']);
            disp('Printed');
        catch,
        end;
    end;
end;