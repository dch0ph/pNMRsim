basename='pmlgexactrfs1';
eval(['load ' basename '_rawamps; amps=' basename '_rawamps;']);
eval(['load ' basename '_rawphases; phases=' basename '_rawphases;']);

subplot(2,1,1);
plot(amps);
subplot(2,1,2);
plot(mod(phases,360));
