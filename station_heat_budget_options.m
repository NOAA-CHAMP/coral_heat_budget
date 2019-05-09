function opts = station_heat_budget_options(stn_or_stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,substitute_field_names)
%function opts = station_heat_budget_options(stn_or_stnm,RAPFX,KMPFX,ISPFX,TIDEPFX,WAVEPFX,substitute_field_names)
%
% Set initial site-specific options for ocean heat budget. Usually called
% from, e.g., STATION_HEAT_BUDGET (see).
%
% Last Saved Time-stamp: <Sun 2014-09-28 18:28:15 Eastern Daylight Time gramer>

  stn = get_station_from_station_name(stn_or_stnm);
  stnm = lower(stn.station_name);
  clear stn_or_stnm;

  station_heat_budget_field_names;

  opts = [];
  % Allow user to preset some options before calling us
  if ( isfield(stn,'opts') )
    opts = stn.opts;
  end;

  switch (stnm),

   case {'lkwf1'},
    % Must be optimized without advection (insufficient data from weekly 1km!)
    opts.default_salinity = get_opt(opts,'default_salinity',36);
    opts.grid_interp_method = get_opt(opts,'grid_interp_method',{'nanmean',4,2});

    opts.override_qlh_adj = get_opt(opts,'override_qlh_adj',0.90);
    opts.do_warm_layer = get_opt(opts,'do_warm_layer',false);
    opts.kd = get_opt(opts,'kd',[0.050,0.400, 67]);
    opts.advection_factor = get_opt(opts,'advection_factor',[0.00,1.00, 45]);
    opts.tidal_advection = get_opt(stn.opts,'tidal_advection',true);
    opts.K_theta = get_opt(opts,'K_theta',[0,20, 45]);
    opts.hc_scaling = get_opt(opts,'hc_scaling','SU');
    opts.hc_warming_factor = get_opt(opts,'hc_warming_factor',0.66);

   case {'fwyf1'},
    opts.default_salinity = get_opt(opts,'default_salinity',35.5);
    % % Results of running OPTIM_Q0 on ERAI met data
    % % *NOTE*: With *in situ* met data, Warm Layer and different/higher Kd are better!
    % opts.kd = get_opt(opts,'kd',[0.10,0.30,0]);
    % opts.do_warm_layer = get_opt(opts,'do_warm_layer',false);
    % % opts.b_convective_coefficient = get_opt(opts,'b_convective_coefficient',8.0e-4);
    % % opts.hc_scaling = get_opt(opts,'hc_scaling','US');
    opts.do_warm_layer = get_opt(opts,'do_warm_layer',true);
    opts.kd = get_opt(opts,'kd',[0.025,0.225,102]);
    opts.advection_factor = get_opt(opts,'advection_factor',[0.00,0.50,102]);
    opts.tidal_advection = get_opt(stn.opts,'tidal_advection',true);
    opts.K_theta = get_opt(opts,'K_theta',[0,20, 45]);
    opts.hc_scaling = get_opt(opts,'hc_scaling','SU');
    opts.hc_warming_factor = get_opt(opts,'hc_warming_factor',0.66);

   case {'mlrf1'},
    opts.default_salinity = get_opt(opts,'default_salinity',35.5);
    % % Results of running OPTIM_Q0 on ERAI met data
    % % *NOTE*: With *in situ* met data, Warm Layer and different/higher Kd are better!
    % % % Good (in situ)
    % % opts.kd = get_opt(opts,'kd',[0.10,0.30,0]);
    % % % Good (ERAI)
    % % opts.kd = get_opt(opts,'kd',[0.15,0.30,0]);
    % % % Better (ERAI)
    % % opts.kd = get_opt(opts,'kd',[0.05,0.35,45]);
    % % Best (ERAI)
    % opts.kd = get_opt(opts,'kd',[0.045,0.375,45]);
    % opts.do_warm_layer = get_opt(opts,'do_warm_layer',false);
    % % opts.b_convective_coefficient = get_opt(opts,'b_convective_coefficient',8.0e-4);
    % opts.K_theta = get_opt(opts,'K_theta',0);
    opts.do_warm_layer = get_opt(opts,'do_warm_layer',true);
    % % opts.kd = get_opt(opts,'kd',[.038,.250, 67]);
    % opts.kd = get_opt(opts,'kd',[.035,.250, 70]);
    opts.kd = get_opt(opts,'kd',[.035,.250, 69]);  % For consistency with SMKF1
    opts.advection_factor = get_opt(opts,'advection_factor',[0.00,1.00, 45]);
    opts.tidal_advection = get_opt(stn.opts,'tidal_advection',true);
    opts.K_theta = get_opt(opts,'K_theta',[0,20, 45]);
    opts.hc_scaling = get_opt(opts,'hc_scaling','SU');
    opts.hc_warming_factor = get_opt(opts,'hc_warming_factor',0.66);

   case {'lonf1'},
    opts.default_salinity = get_opt(opts,'default_salinity',35.7);
    % %%%%DEBUG???
    % % opts.grid_interp_method = get_opt(opts,'grid_interp_method','triangular,linear');
    % % Results of running OPTIM_Q0 on ERAI met data
    % % *NOTE*: With *in situ* met data, Warm Layer and different/higher Kd are better!
    % % opts.kd = get_opt(opts,'kd',[0.30,0.50,274]);
    % opts.kd = get_opt(opts,'kd',[0.50,0.70,0]);
    % opts.do_warm_layer = get_opt(opts,'do_warm_layer',false);
    % % opts.b_convective_coefficient = get_opt(opts,'b_convective_coefficient',8.0e-4);
    opts.do_warm_layer = get_opt(opts,'do_warm_layer',false);
    % opts.kd = get_opt(opts,'kd',[0.500,1.250, 45]);
    opts.kd = get_opt(opts,'kd',[.475,1.275, 45]);
    opts.advection_factor = get_opt(opts,'advection_factor',[0.00,1.00, 45]);
    opts.tidal_advection = get_opt(stn.opts,'tidal_advection',false);
    opts.K_theta = get_opt(opts,'K_theta',[0,20, 45]);
    opts.hc_scaling = get_opt(opts,'hc_scaling','SU');
    opts.hc_warming_factor = get_opt(opts,'hc_warming_factor',0.66);

   case {'smkf1'},
    opts.default_salinity = get_opt(opts,'default_salinity',35.5);
    % % Results of running OPTIM_Q0 on ERAI *AND* in situ met data
    % % *NOTE*: Earlier peak jd=91 improves seasonal phase problem a lot!
    % %opts.kd = get_opt(opts,'kd',[0.10,0.30,91]);
    % opts.kd = get_opt(opts,'kd',[0.10,0.30,274]);
    % opts.do_warm_layer = get_opt(opts,'do_warm_layer',true);
    % % opts.b_convective_coefficient = get_opt(opts,'b_convective_coefficient',8.0e-4);
    opts.do_warm_layer = get_opt(opts,'do_warm_layer',true);
    % 0.100,0.500, 33
    opts.kd = get_opt(opts,'kd',[0.066,0.450, 69]);
    opts.advection_factor = get_opt(opts,'advection_factor',[0.00,1.00, 45]);
    opts.tidal_advection = get_opt(stn.opts,'tidal_advection',true);
    opts.K_theta = get_opt(opts,'K_theta',[0,20, 45]);
    opts.hc_scaling = get_opt(opts,'hc_scaling','UNSTEADY_SEASONAL');
    opts.hc_warming_factor = get_opt(opts,'hc_warming_factor',0.66);

   case {'looe1'},
    opts.default_salinity = get_opt(opts,'default_salinity',36);
    % Results of running OPTIM_Q0 on ERAI met data
    % opts.b_convective_coefficient = get_opt(opts,'b_convective_coefficient',8.0e-4);
    switch ( lower(sfld) ),
     case {'microcat_seatemp','mc_seatemp'},
      % opts.kd = get_opt(opts,'kd',0.3);
      % opts.do_warm_layer = get_opt(opts,'do_warm_layer',true);
      % opts.K_theta = get_opt(opts,'K_theta',10);
      opts.kd = get_opt(opts,'kd',[0.050,0.150, 45]);
      opts.do_warm_layer = get_opt(opts,'do_warm_layer',true);
      opts.advection_factor = get_opt(opts,'advection_factor',[0.00,0.25, 45]);
      opts.tidal_advection = get_opt(stn.opts,'tidal_advection',true);
      opts.K_theta = get_opt(opts,'K_theta',[0,20, 45]);
      opts.hc_scaling = get_opt(opts,'hc_scaling','UNSTEADY_SEASONAL');
      opts.hc_warming_factor = get_opt(opts,'hc_warming_factor',0.66);
     case {'adcp_seatemp','ad_seatemp'},
      opts.keep_bad_dates = get_opt(opts,'keep_bad_dates',true);
      opts.kd = get_opt(opts,'kd',[0.025,0.200, 80]);
      opts.do_warm_layer = get_opt(opts,'do_warm_layer',true);
      opts.advection_factor = get_opt(opts,'advection_factor',[0.00,0.25, 45]);
      opts.tidal_advection = get_opt(stn.opts,'tidal_advection',true);
      opts.K_theta = get_opt(opts,'K_theta',[0,20, 45]);
      opts.hc_scaling = get_opt(opts,'hc_scaling','STEADY_SEASONAL');
      opts.hc_warming_factor = get_opt(opts,'hc_warming_factor',0.66);
     otherwise,
      error('Do not yet handle LOOE with sea temperature "%s"',sfld);
    end;

   case {'sanf1'},
    opts.default_salinity = get_opt(opts,'default_salinity',36);
    opts.grid_interp_method = get_opt(opts,'grid_interp_method','nearest');

    % % Results of running OPTIM_Q0 on ERAI *AND* in situ met data
    % % *NOTE*: Earlier peak jd=91 improves seasonal phase problem slightly
    % % *NOTE*: But no matter what, all results at this site are *TERRIBLE*!
    % %opts.kd = get_opt(opts,'kd',[0.10,0.30,91]);
    % opts.kd = get_opt(opts,'kd',[0.10,0.30,274]);
    % opts.do_warm_layer = get_opt(opts,'do_warm_layer',false);
    % opts.b_convective_coefficient = get_opt(opts,'b_convective_coefficient',8.0e-4);
    opts.do_warm_layer = get_opt(opts,'do_warm_layer',true);
    opts.kd = get_opt(opts,'kd',[0.015,0.150, 67]);
    opts.advection_factor = get_opt(opts,'advection_factor',0);
    opts.tidal_advection = get_opt(stn.opts,'tidal_advection',true);
    opts.K_theta = get_opt(opts,'K_theta',[0,20, 45]);
    opts.hc_scaling = get_opt(opts,'hc_scaling','UNSTEADY_SEASONAL');
    opts.hc_warming_factor = get_opt(opts,'hc_warming_factor',0.66);

   case {'dryf1'},
    opts.default_salinity = get_opt(opts,'default_salinity',36);
    opts.grid_interp_method = get_opt(opts,'grid_interp_method',{'nanmean',2,2});

    % % Results of running OPTIM_Q0 on ERAI *AND* in situ met data
    % % *NOTE*: Earlier peak jd=91 improves seasonal phase problem slightly
    % % *NOTE*: But no matter what, all results at this site are *TERRIBLE*!
    % %opts.kd = get_opt(opts,'kd',[0.10,0.30,91]);
    % opts.kd = get_opt(opts,'kd',[0.10,0.30,274]);
    % opts.do_warm_layer = get_opt(opts,'do_warm_layer',false);
    % opts.b_convective_coefficient = get_opt(opts,'b_convective_coefficient',8.0e-4);
    opts.override_qlh_adj = get_opt(opts,'override_qlh_adj',0.90);
    opts.do_warm_layer = get_opt(opts,'do_warm_layer',false);
    opts.kd = get_opt(opts,'kd',[0.150,0.500,354]);
    opts.advection_factor = get_opt(opts,'advection_factor',[0.00,1.00, 45]);
    opts.tidal_advection = get_opt(stn.opts,'tidal_advection',true);
    opts.K_theta = get_opt(opts,'K_theta',[0,20, 45]);
    opts.hc_scaling = get_opt(opts,'hc_scaling','SU');
    opts.hc_warming_factor = get_opt(opts,'hc_warming_factor',0.66);

  end;


  %%% GENERAL DEFAULTS (if not already overridden for a specific station above)

  % % opts.kd = get_opt(opts,'kd',[0.10,0.25,274]);
  % opts.kd = get_opt(opts,'kd',[0.10,0.30,274]);
  % % opts.kd = get_opt(opts,'kd',[0.10,0.40,274]);
  % % % Contrary to evidence
  % % opts.kd = get_opt(opts,'kd',[0.10,0.40,90]);

  % "DEFAULT": 0.20
  opts.kd = get_opt(opts,'kd',0.20);
  % opts.kd = get_opt(opts,'kd',0.40);

  % % Try NOT applying warm layer adjustment
  % opts.do_warm_layer = get_opt(opts,'do_warm_layer',false);

  % % %%%% ??? DEBUG - try using minimum bottom reflectivity
  % % opts.sand_fraction = get_opt(opts,'sand_fraction',0.00);
  % % %%%% ??? DEBUG - try "total" attenuation
  % % opts.absorption_factor_gamma = get_opt(opts,'absorption_factor_gamma',1.0);
  % %%%% ??? DEBUG - try turning off bottom flux
  % opts.bottom_reflectance = get_opt(opts,'bottom_reflectance',1.0);

  % opts.b_convective_coefficient = get_opt(opts,'b_convective_coefficient',2e-4);

  % opts.add_alongshore_advection = get_opt(opts,'add_alongshore_advection',false);
  % opts.gradient_climatology = get_opt(opts,'gradient_climatology',[-2e-4, +2e-4, 355]);
  %%%% ??? DEBUG
  % opts.gradient_climatology = get_opt(opts,'gradient_climatology',[0, 0, 355]);


  % % Based on Laplacians of 1km satellite SST and model surface T
  % opts.laplacian_climatology = get_opt(opts,'laplacian_climatology',[-1.3e-8, +2.2e-8, 183]);
  % % Center seasonal Laplacian at 0? trying to eliminate annual bias!
  % opts.laplacian_climatology = get_opt(opts,'laplacian_climatology',[-2.2e-8, +2.2e-8, 183]);
  %%%% ??? DEBUG
  % opts.laplacian_climatology = get_opt(opts,'laplacian_climatology',[0, 0, 183]);


  % opts.K_theta = get_opt(opts,'K_theta',[model_K_theta,model_K_theta,90]);
  % % opts.K_theta = get_opt(opts,'K_theta',[2,20,1]);
  % % opts.K_theta = get_opt(opts,'K_theta',[20,20,1]);
  % % opts.K_theta.func = @(ts)(struct('date',ts.date,'data',repmat(20,size(ts.data))));
  % % opts.K_theta.func = @(ts)(struct('date',ts.date,'data',2.0+(0.5.*ts.data)));
  % % opts.K_theta.arg = Wlpfld;


  % % Try ignoring km-scale (hydrodynamic model) advection
  % opts.km_scale_advection = get_opt(opts,'km_scale_advection',false);

  % % Try ignoring km-scale sea temperature gradient
  % opts.gradient_climatology = get_opt(opts,'gradient_climatology',0);

  % % Try ignoring km-scale eddy thermal diffusion
  % opts.laplacian_climatology = get_opt(opts,'laplacian_climatology',0);



  % opts.hc_scaling = get_opt(opts,'hc_scaling','Adaptive');
  % opts.hc_scaling = get_opt(opts,'hc_scaling','Seasonal');
  % opts.hc_scaling = get_opt(opts,'hc_scaling','SS');
  % opts.hc_scaling = get_opt(opts,'hc_scaling','US');
  opts.hc_scaling = get_opt(opts,'hc_scaling','SU');
  % opts.hc_scaling = get_opt(opts,'hc_scaling','UU');

  % opts.hc_warming_factor = get_opt(opts,'hc_warming_factor',0.9);
  % opts.hc_cooling_factor = get_opt(opts,'hc_cooling_factor',1.2);
  % opts.hc_warming_factor = get_opt(opts,'hc_warming_factor',1.3);
  % opts.hc_cooling_factor = get_opt(opts,'hc_cooling_factor',2.0);

  disp(opts);

return;
