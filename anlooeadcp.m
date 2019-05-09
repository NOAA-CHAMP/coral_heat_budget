1;

error('Old script');

looe = load('C:\Documents and Settings\gramer\My Documents\RSMAS\Coastal\thesis\data\SFP\Looe_Key_ADCP_21JUN2008_thru_25JAN2009.mat');

dts = datenum(looe.SerYear,looe.SerMon,looe.SerDay,looe.SerHour,looe.SerMin,looe.SerSec);

bn=2;
figure;
plot(dts,[looe.SerEmmpersec(:,bn),looe.SerNmmpersec(:,bn),looe.SerVmmpersec(:,bn)]/1000);
maxigraph; datetick3;
titlename(['Bin ' num2str(bn)]);
legend('u','v','w', 'Location','Best');

bn=34;
figure;
plot(dts,[looe.SerEmmpersec(:,bn),looe.SerNmmpersec(:,bn),looe.SerVmmpersec(:,bn)]/1000);
maxigraph; datetick3;
titlename(['Bin ' num2str(bn)]);
legend('u','v','w', 'Location','Best');


bn1=34;bn2=2; figure; plot(dts,[looe.SerEmmpersec(:,34)-looe.SerEmmpersec(:,2),looe.SerNmmpersec(:,34)-looe.SerNmmpersec(:,2)])/1000); maxigraph; datetick3; titlename(['\Delta V (' num2str(bn1) ',' num2str(bn2) ')']); legend('u','v', 'Location','Best');
bn1=34;bn2=2; figure; plot(dts,[looe.SerEmmpersec(:,bn1)-looe.SerEmmpersec(:,bn2),looe.SerNmmpersec(:,bn1)-looe.SerNmmpersec(:,bn2])/1000); maxigraph; datetick3; titlename(['\Delta V (' num2str(bn1) ',' num2str(bn2) ')']); legend('u','v', 'Location','Best');
bn1=34; bn2=2; figure; plot(dts,([looe.SerEmmpersec(:,bn1)-looe.SerEmmpersec(:,bn2),looe.SerNmmpersec(:,bn1)-looe.SerNmmpersec(:,bn2)])/1000); maxigraph; datetick3; titlename(['\Delta V (' num2str(bn1) ',' num2str(bn2) ')']); legend('u','v', 'Location','Best');

%%%% looe1 = plot_ngdc_bathy_station('looe1');
bn1=34;
bn2=2;
figure;
plot(dts,abs([looe.SerEmmpersec(:,bn1)-looe.SerEmmpersec(:,bn2),looe.SerNmmpersec(:,bn1)-looe.SerNmmpersec(:,bn2)])/1000);
maxigraph; datetick3;
titlename(['|V' num2str(bn1) ' - V' num2str(bn2) '|']); legend('u','v', 'Location','Best');

bn1=34;
bn2=2;
figure;
plot(dts,([looe.SerDir10thDeg(:,bn1)-looe.SerDir10thDeg(:,bn2)])/10);
maxigraph; datetick3;
titlename(['Dir (' num2str(bn1) '-' num2str(bn2) ')']);
