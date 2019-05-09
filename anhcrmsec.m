1;

hcrmsec=[];
flds={};
for cfld=fieldnames(s)';
  fld=cfld{:};
  if ( isfield(s.(fld),'hcrmsec') )
    disp(sprintf('% 8.2f\t%s',s.(fld).hcrmsec,fld));
    hcrmsec(end+1)=s.(fld).hcrmsec;
    flds{end+1}=fld;
  end;
end

fmg;
plot(hcrmsec,1:numel(hcrmsec),'.-');
set(gca,'ytick',1:numel(hcrmsec),'yticklabel',flds,'xscale','log');
xlim([0.1,200]);
titlename([upper(s.(fld).station_name),' Sensitivity Analysis']);
print('-dtiff',fullfile(get_thesis_path('../figs'),[lower(s.(fld).station_name),'-anhcrmsec.tif']));
