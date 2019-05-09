function plot_adcp_bins(stn,bn1,bn2,plotflgs)
%function plot_adcp_bins(stn,[bn1],[bn2],[plotflgs])
%
% Last Saved Time-stamp: <Tue 2010-11-16 14:04:26  Lew.Gramer>

  if ( ~exist('bn1','var') || isempty(bn1) )
    bn1 = 25;
  end;
  if ( ~exist('bn2','var') || isempty(bn2) )
    bn2 = 2;
  end;
  if ( ~exist('plotflgs','var') || isempty(plotflgs) )
    plotflgs = [0 0 1 0 0];
  end;
  if ( isscalar(plotflgs) )
    plotflgs = [plotflgs plotflgs plotflgs plotflgs plotflgs];
  end;

  if plotflgs(1)
    figure;
    plot(stn.adcp_u.date,[stn.adcp_u.prof(:,bn2),stn.adcp_v.prof(:,bn2),...
                        stn.adcp_w.prof(:,bn2)].*1e3);
    maxigraph; datetick3;
    titlename(['Bin ' num2str(bn2)]);
    legend('u','v','w', 'Location','Best');
    ylabel('[ mm / s ]');
    hold on;
  end;

  if plotflgs(2)
    figure;
    plot(stn.adcp_u.date,[stn.adcp_u.prof(:,bn1),stn.adcp_v.prof(:,bn1),...
                        stn.adcp_w.prof(:,bn1)].*1e3);
    maxigraph; datetick3;
    titlename(['Bin ' num2str(bn1)]);
    legend('u','v','w', 'Location','Best');
    ylabel('[ mm / s ]');
    hold on;
  end;

  if plotflgs(3)
    figure;
    plot(stn.adcp_u.date,([stn.adcp_u.prof(:,bn1)-stn.adcp_u.prof(:,bn2),...
                        stn.adcp_v.prof(:,bn1)-stn.adcp_v.prof(:,bn2)]).*1e3);
    maxigraph; datetick3;
    titlename(['\Delta V (' num2str(bn1) ',' num2str(bn2) ')']);
    legend('u','v', 'Location','Best');
    ylabel('[ mm / s ]');
    hold on;
  end;

  if plotflgs(4)
    figure;
    plot(stn.adcp_u.date,abs([stn.adcp_u.prof(:,bn1)-stn.adcp_u.prof(:,bn2),...
                        stn.adcp_v.prof(:,bn1)-stn.adcp_v.prof(:,bn2)]).*1e3);
    maxigraph; datetick3;
    titlename(['|V' num2str(bn1) ' - V' num2str(bn2) '|']); legend('u','v', 'Location','Best');
    ylabel('[ mm / s ]');
    hold on;
  end;

  if plotflgs(5)
    figure;
    plot(stn.adcp_u.date,([stn.adcp_dir.prof(:,bn1)-stn.adcp_dir.prof(:,bn2)]));
    maxigraph; datetick3;
    titlename(['Dir (' num2str(bn1) '-' num2str(bn2) ')']);
    hold on;
  end;

return;
