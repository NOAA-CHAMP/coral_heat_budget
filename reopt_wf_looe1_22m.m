1;

diary(fullfile(get_thesis_path('../src'),'reopt_wf_looe1_22m.log'));
more off
timenow

%{
disp('**********************************************************************');
disp('CONTROL');
disp('**********************************************************************');
stn = optimize_station_heat_budget('looe1','erai','avhrr_weekly','ndbc','tpxo_tide','erai',{'sfld','ad_seatemp'});
stn=[]; clear stn;

    {'_24_h_lp','UNSTEADY_SEASONAL',0.66,{[0.050,0.150, 45],},{[0.00,0.25, 45]},{[0,20, 45]},true},... %err  5.2,1.2
    {'_24_h_lp','STEADY_SEASONAL',  0.66,{[0.050,0.150, 45],},{[0.00,0.25, 45]},{[0,20, 45]},true},... %err  2.1,1.3
    {'_24_h_lp','UNSTEADY_SEASONAL',0.66,{[0.050,0.150, 69],},{[0.00,0.25, 45]},{[0,20, 45]},true},... %err  5.7,1.2
    {'_24_h_lp','STEADY_SEASONAL',  0.66,{[0.050,0.150, 69],},{[0.00,0.25, 45]},{[0,20, 45]},true},... %err  ?.2,1.2
    {'_24_h_lp','UNSTEADY_SEASONAL',0.66,{[0.050,0.200, 45],},{[0.00,0.25, 45]},{[0,20, 45]},true},... %err  ?.2,1.2
    {'_24_h_lp','STEADY_SEASONAL',  0.66,{[0.050,0.200, 45],},{[0.00,0.25, 45]},{[0,20, 45]},true},... %err  ?.2,1.2
    {'_24_h_lp','UNSTEADY_SEASONAL',0.66,{[0.050,0.200, 69],},{[0.00,0.25, 45]},{[0,20, 45]},true},... %err  ?.2,1.2
    {'_24_h_lp','STEADY_SEASONAL',  0.66,{[0.050,0.200, 69],},{[0.00,0.25, 45]},{[0,20, 45]},true},... %err  ?.2,1.2

%%%% BETTER
    {'_24_h_lp','STEADY_SEASONAL',  0.66,{[0.050,0.150,113],},{[0.00,0.25, 45]},{[0,20, 45]},true},... %err  3.4,1.7
    {'_24_h_lp','STEADY_SEASONAL',  0.66,{[0.025,0.250,113],},{[0.00,0.25, 45]},{[0,20, 45]},true},... %err  2.4,1.7

%%%% NEARING OPTIMAL
    {'_24_h_lp','STEADY_SEASONAL',  0.66,{[0.025,0.150,126],},{[0.00,0.25, 45]},{[0,20, 45]},true},... %err  2.4,1.4
    {'_24_h_lp','STEADY_SEASONAL',  0.66,{[0.025,0.200,126],},{[0.00,0.25, 45]},{[0,20, 45]},true},... %err  2.0,1.7
    {'_24_h_lp','STEADY_SEASONAL',  0.66,{[0.025,0.150,101],},{[0.00,0.25, 45]},{[0,20, 45]},true},... %err  2.9,1.1
    {'_24_h_lp','STEADY_SEASONAL',  0.66,{[0.025,0.200,101],},{[0.00,0.25, 45]},{[0,20, 45]},true},... %err  1.9,1.3
    {'_24_h_lp','STEADY_SEASONAL',  0.66,{[0.025,0.150, 92],},{[0.00,0.25, 45]},{[0,20, 45]},true},... %err  3.3,1.0
%}

%%%% CONTROL: SU,wf=1,[0.01,0.40,113],[0,0.25,45],[0,20,45],false: err 2.1, 1.9
for cop={ ...
    %%%{'_24_h_lp','STEADY_SEASONAL',  0.66,{[0.050,0.150, 45],},{[0.00,0.25, 45]},{[0,20, 45]},true},... %err  2.1,1.3
    %%%{'_24_h_lp','STEADY_SEASONAL',  0.66,{[0.025,0.150,113],},{[0.00,0.25, 45]},{[0,20, 45]},true},... %err  2.6,1.2
    %%%{'_24_h_lp','STEADY_SEASONAL',  0.66,{[0.025,0.200,113],},{[0.00,0.25, 45]},{[0,20, 45]},true},... %err  1.9,1.5
    %%{'_24_h_lp','STEADY_SEASONAL',  0.66,{[0.025,0.200, 92],},{[0.00,0.25, 45]},{[0,20, 45]},true},... %err  1.9,1.1
    %%{'_24_h_lp','STEADY_SEASONAL',  0.66,{[0.025,0.150, 80],},{[0.00,0.25, 45]},{[0,20, 45]},true},... %err  3.8,0.9
    %%{'_24_h_lp','STEADY_SEASONAL',  0.66,{[0.025,0.200, 80],},{[0.00,0.25, 45]},{[0,20, 45]},true},... %err  2.2,1.0
    %%{'_24_h_lp','STEADY_SEASONAL',  0.66,{[0.025,0.175, 80],},{[0.00,0.25, 45]},{[0,20, 45]},true},... %err  2.8,0.9
    %%{'_24_h_lp','STEADY_SEASONAL',  0.66,{[0.025,0.175, 86],},{[0.00,0.25, 45]},{[0,20, 45]},true},... %err  2.5,1.0
    %{'_24_h_lp','STEADY_SEASONAL',  0.66,{[0.025,0.200, 92],},{[0.00,0.25, 45]},{[0,20, 45]},true},... %err  1.9,1.1
    %{'_24_h_lp','STEADY_SEASONAL',  0.66,{[0.025,0.190, 80],},{[0.00,0.25, 45]},{[0,20, 45]},true},... %err  2.4,0.9
    %{'_24_h_lp','STEADY_SEASONAL',  0.66,{[0.025,0.190, 82],},{[0.00,0.25, 45]},{[0,20, 45]},true},... %err  2.3,1.0
    %{'_24_h_lp','STEADY_SEASONAL',  0.66,{[0.025,0.210, 80],},{[0.00,0.25, 45]},{[0,20, 45]},true},... %err  2.0,1.0
    %{'_24_h_lp','STEADY_SEASONAL',  0.66,{[0.025,0.210, 82],},{[0.00,0.25, 45]},{[0,20, 45]},true},... %err  2.0,1.0
    {'_24_h_lp','STEADY_SEASONAL',  0.66,{[0.025,0.200, 80],},{[0.00,0.25, 45]},{[0,20, 45]},true},... %err  2.2,1.0
        };

    op = cop{:};
    lp = op{1};
    sc = op{2};
    wf = op{3};
    kd = op(4);
    ad = op(5);
    kt = op(6);
    kb = op{7};
    disp('**********************************************************************');
    disp(lp);
    disp({sc,wf,kd});
    disp(kd{:});
    disp(ad{:});
    disp(kt{:});
    disp(kb);
    disp('**********************************************************************');

    stn = optimize_station_heat_budget('looe1','erai','avhrr_weekly','ndbc','tpxo_tide','erai',{'sfld','ad_seatemp','Q0_LOWPASS',lp},struct('hc_scaling',sc,'hc_warming_factor',wf,'override_light_kds',kd,'override_advection_factors',ad,'override_k_thetas',kt,'keep_bad_dates',kb));
    stn=[]; clear stn;
    pause(1);

end;

timenow
more on
diary off