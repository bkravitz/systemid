function [] = plot_time_freq(input)


figure;
subplot(1,2,1);
plot(input);
title('Time Domain','FontSize',15);
ylabel('Temperature Perturbation (K)','FontSize',15);
set(gca,'FontSize',15,'XTick',[365*1 365*2 365*3 365*4 365*5],'XTickLabel',{'1yr','2yr','3yr','4yr','5yr'});
xlabel('Time','FontSize',15);
% xlim([1 7300]);
% ylim([-1 1]);

subplot(1,2,2)
[Pxx,F] = pwelch(input);
semilogx(F/2/pi, 20*log10(Pxx));
set(gca, 'XTick',sort([1/2 1/7 1/30 1/91 1/182 1/365 1/365/2 1/365/5]),...
    'XTickLabel',{'5yr','2yr','1yr','6mo','3mo','1mo','1wk','2d'});
ylabel('dB','FontSize',15);
xlabel('Timescale Corresponding to (Angular) Frequency','FontSize',15);
set(gca,'FontSize',15);
% xlim([1/365/10 1/2]);
title('Frequency Domain','FontSize',15);

end