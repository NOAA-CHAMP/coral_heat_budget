function v = get_nc_insol(fname)

  vars = { ...
      'V_delta_lat', ...
      'V_delta_lon', ...
      ...
      'V_lat_cell', ...
      'V_lon_cell', ...
      ...
      'V_frac_total_cld', ...
      'V_lwdsfc', ...
      'V_lwusfc', ...
      'V_olr', ...
      'V_swd_sfc', ...
      'V_swu_sfc', ...
      'V_vis_down_sfc', ...
         };

  v.dlat = ncload(fname, 'V_delta_lat');
  lats = ncload(fname, 'V_lat_cell');
  v.dlon = ncload(fname, 'V_delta_lon');
  lons = ncload(fname, 'V_lon_cell');

  [LONS,LATS] = meshgrid(lons, lats);
  [STN_LONS,STN_LATS] = meshgrid(lons, lats);

  idx = find(stn_lons
  v.lats = 

return;


% f =
% 
%           NetCDF_File: '/phoddat/share/gramer/nesdis/satpar/foo.nc'
%           nDimensions: 7
%            nVariables: 12
%     nGlobalAttributes: 5
%       RecordDimension: 'time'
%              nRecords: 0
%            Permission: 'write'
%            DefineMode: 'data'
%              FillMode: 'fill'
%            MaxNameLen: 0



% nc = netcdf('foo.nc', 'noclobber');
% if isempty(nc), return, end

% %% Global attributes:

% nc.Description = ncchar(''Test file.'');
% nc.Author = ncchar(''Dr. Charles R. Denham, ZYDECO.'');
% nc.Created = ncchar(''27-Jan-2010 09:01:21'');
% nc.a.dotted.name = ncchar(''a dotted global attribute name'');
% nc.a-dashed-name = ncchar(''a dashed global attribute name'');

% %% Dimensions:

% nc('time') = 0; %% (record dimension)
% nc('lat') = 10;
% nc('lon') = 5;
% nc('elapsed_time') = 100;
% nc('horse_number') = 5;
% nc('a.dotted.name') = 22;
% nc('a-dashed-name') = 22;

% %% Variables and attributes:

% nc{'time'} = nclong('time'); %% 0 elements.
% nc{'time'}.units = ncchar(''seconds'');

% nc{'lat'} = nclong('lat'); %% 10 elements.
% nc{'lat'}.FillValue_ = nclong(-999);
% nc{'lat'}.scale_factor = nclong(2);
% nc{'lat'}.add_offset = nclong(100);
% nc{'lat'}.units = ncchar(''degrees_north'');

% nc{'lon'} = nclong('lon'); %% 5 elements.
% nc{'lon'}.units = ncchar(''degrees_east'');

% nc{'elapsed_time'} = ncdouble('elapsed_time'); %% 100 elements.
% nc{'elapsed_time'}.units = ncchar(''fortnights'');

% nc{'horse_number'} = nclong('horse_number'); %% 5 elements.

% nc{'speed'} = ncdouble('elapsed_time', 'horse_number'); %% 500 elements.
% nc{'speed'}.units = ncchar(''furlongs/fortnight'');

% nc{'z'} = nclong('time', 'lat', 'lon'); %% 0 elements.
% nc{'z'}.units = ncchar(''meters'');
% nc{'z'}.valid_range = ncfloat([0 5000]);

% nc{'t'} = nclong('time', 'lat', 'lon'); %% 0 elements.

% nc{'p'} = nclong('time', 'lat', 'lon'); %% 0 elements.
% nc{'p'}.FillValue_ = nclong(-1);

% nc{'rh'} = nclong('time', 'lat', 'lon'); %% 0 elements.
% nc{'rh'}.FillValue_ = nclong(-1);

% nc{'a.dotted.name'} = ncchar('a.dotted.name'); %% 22 elements.
% nc{'a.dotted.name'}.a.dotted.name = ncchar(''a dotted variable attribute name'');

% nc{'a-dashed-name'} = ncchar('a-dashed-name'); %% 22 elements.
% nc{'a-dashed-name'}.a-dashed-name = ncchar(''a dashed variable attribute name'');

% endef(nc)
% close(nc)
