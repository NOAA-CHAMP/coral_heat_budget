function [baseurl,sat,yr,mo,dy,hr,mn,region,rest,fmt] = parseusfurl(url)
%function [baseurl,sat,yr,mo,dy,hr,mn,region,rest,fmt] = parseusfurl(url)
%
% Sample: http://.../n12.20060601.0948.florida.true.png

  %[baseurl,fname,fmt,ig] = fileparts(url);
  [baseurl,fname,fmt] = fileparts(url);

  Cstr = textscan(fname, '%[^.].%04d%02d%02d.%02d%02d.%[^.].%s');

  sat = Cstr{1}{:};

  yr = double(Cstr{2});
  mo = double(Cstr{3});
  dy = double(Cstr{4});
  hr = double(Cstr{5});
  mn = double(Cstr{6});

  region = Cstr{7}{:};
  rest = Cstr{8}{:};

return;
