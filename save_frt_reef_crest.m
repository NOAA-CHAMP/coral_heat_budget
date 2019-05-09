1;

cresthi=-11;
crestlo=-60;

fh=fmg; [ig,ig,cshi,ch] = map_freef([-81.9,-80.0,24.2,26.0],cresthi); close(fh);
fh=fmg; [ig,ig,cslo,ch] = map_freef([-81.9,-80.0,24.2,26.0],crestlo); close(fh);

bkix=find(cshi(1,:)==cresthi); cs = [cslo(:,end:-1:2) , cshi(:,bkix(1)+1:1:bkix(2)-1)];

crestlon=cs(1,:);
crestlat=cs(2,:);

save('frt_reef_crest.mat','crestlon','crestlat','cresthi','crestlo');
