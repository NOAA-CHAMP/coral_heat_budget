function stn = calc_gom_hycom_terms(stn,dtix,doPlot,fld)
%function stn = calc_gom_hycom_terms(stn,dtix,doPlot,fld)

  if ( ~exist('dtix','var') || isempty(dtix) )
    dtix = 3;
  end;
  if ( ~exist('doPlot','var') || isempty(doPlot) )
    doPlot = true;
  end;
  if ( ~exist('fld','var') || isempty(fld) )
    fld = 'gom_hycom_seatemp_field';
  end;

  dts = stn.(fld).date;
  lon = stn.(fld).lon;
  lat = stn.(fld).lat;
  fld = stn.(fld).field;


  %DEBUG:  stn.(fld)=rmfield(stn.(fld),{'gradient_x','gradient_y','gradient_t','laplacian'});

  if ( all(isfield(stn.(fld),{'gradient_x','gradient_y','gradient_t'})) )
    gx = stn.(fld).gradient_x;
    gy = stn.(fld).gradient_y;
    dTdt = stn.(fld).gradient_t;
  else
    % Use MKS units for gradient and Laplacian
    dt = 1;
    dx = sw_dist(lat([1 1]),lon([1 2]),'km')*1e3;
    dy = sw_dist(lat([1 2]),lon([1 1]),'km')*1e3;
    rotfld = permute(fld,[2 3 1]);

    [gx,gy,dTdt] = gradient(rotfld,dx,dy,dt);
    gx = permute(gx,[3 1 2]);
    gy = permute(gy,[3 1 2]);
    dTdt = permute(dTdt,[3 1 2]);
    stn.(fld).gradient_x = gx;
    stn.(fld).gradient_y = gy;
    stn.(fld).gradient_t = dTdt;
  end;

  if ( isfield(stn.(fld),'laplacian') )
    l = stn.(fld).laplacian;
  else
    [nt,nx,ny] = size(fld);
    l = [];
    for ix=1:nt
      l(ix,1:nx,1:ny) = del2(squeeze(fld(ix,:,:)),dx,dy);
    end;
    stn.(fld).laplacian = l;
  end;

  if ( doPlot )
    figure; maxigraph;
    ax1=subplot(2,2,1); contourf(lon,lat,squeeze(fld(dtix,:,:))); title('T');
    ax2=subplot(2,2,2); contourf(lon,lat,squeeze(gx(dtix,:,:))); title('\nablaT^x');
    ax3=subplot(2,2,4); contourf(lon,lat,squeeze(gy(dtix,:,:))); title('\nablaT^y');
    ax4=subplot(2,2,3); contourf(lon,lat,squeeze(l(dtix,:,:))); title('\nabla^2T');
    suptitle(datestr(dts(dtix)));
    colorbar('peer',ax1); colorbar('peer',ax2); colorbar('peer',ax3); colorbar('peer',ax4);
  end;

return;
