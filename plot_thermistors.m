1;

x = csvread('../thermistor.csv',1,1);
bbox = [min(x(:,2)) max(x(:,2)) min(x(:,1)) max(x(:,1))];
map_sofla(bbox,[-10 -30 -200]);
set(gcf,'Units','normalized','OuterPosition',[0 0 1 1]);
hold on;
plot(x(:,2),x(:,1),'p');
text(x(:,2),x(:,1),num2str(x(:,3)));
