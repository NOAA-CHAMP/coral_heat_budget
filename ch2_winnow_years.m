1;

fld='ndbc_sea_t';

warning('OFF','run_length_encode:generic');

ys=1987:2014;
for cst={fwyf1,mlrf1,lonf1,smkf1,sanf1};
  %ys=1987:2014;
  t=cst{:}.(fld);
  disp(cst{:}.station_name);

  % Only use complete days - those with 21+ hourly data points
  uns = floor(t.date);
  [rl,badun] = run_length_encode(uns);
  badun(rl>=22) = [];
  t.data(ismember(uns,badun)) = [];
  t.date(ismember(uns,badun)) = [];

  % Only use complete weeks - those with 5+ full days of data
  uns = get_week(t.date);
  [rl,badun] = run_length_encode(uns);
  badun(rl>=5) = [];
  t.data(ismember(get_week(t.date),badun)) = [];
  t.date(ismember(get_week(t.date),badun)) = [];

  % Only use complete months - those with 24+ full days of data
  uns = get_month(t.date);
  [rl,badun] = run_length_encode(uns);
  badun(rl>=24) = [];
  t.data(ismember(get_month(t.date),badun)) = [];
  t.date(ismember(get_month(t.date),badun)) = [];

  clear uns badun rl

  for yr=ys(:)';
    yrix = find(get_year(t.date)==yr);
    %if (numel(yrix)<(270*24));
    if (numel(yrix)<1);
      %disp(yr);
      ys(ys==yr)=[];
    end;
  end;
  t=[]; cst=[];
  %disp(ys);
end;

warning('ON','run_length_encode:generic');

clear st cst;
disp(ys);
