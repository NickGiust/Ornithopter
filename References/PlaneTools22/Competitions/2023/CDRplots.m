close all; clc;

GMbestList = (20:10:60);

% Score vs. GMbest

scoreList_m3p1_m2p1 = [1.86 1.86 1.97 1.77 1.77];
scoreList_m3p1_m2p2 = [1.89 1.89 2.05 1.85 1.79];
scoreList_m3p1_m2p3 = [2.13 2.13 1.99 1.78 1.78];
scoreList_m3p1_m2p4 = [2.42 2.42 2.32 2.07 2.01];


figure
sz = 100;
scatter(GMbestList, scoreList_m3p1_m2p1, sz, 'filled', 'd')
hold on
scatter(GMbestList, scoreList_m3p1_m2p2, sz, 'filled', 'square')
hold on
scatter(GMbestList, scoreList_m3p1_m2p3, sz, 'filled', '^')
hold on
scatter(GMbestList, scoreList_m3p1_m2p4, sz, 'filled', 'o')
lgd = legend('12x6 6S 4400mAh', '12x8 6S 4400mAh','14x12 4S 6600mAh','14x12 6S 4400mAh');
lgd.FontSize = 18;
lgd.FontName = 'Times New Roman';
xlabel("GM Best Assumption",'FontSize',18,'FontWeight','bold', 'FontName', 'Times New Roman')
xlim([10 70])
ylabel("Score",'FontSize',18,'FontWeight','bold', 'FontName', 'Times New Roman')
ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'off';

%%
% Span vs. GMbest
spanList_m3p1_m2p4 = [2.4 2.4 2 1.8 2.4];
ARList_m3p1_m2p4 = [13 13 11 11 13];
m2payloadFractionList_m3p1_m2p4 = [0.5 0.5 0.5 0.5 0.4];
m3antennaLength_m3p1_m2p4 = [1 1 1 1 1];

figure
subplot(2,1,1)
sz = 80;
scatter(GMbestList, spanList_m3p1_m2p4, sz, 'filled', 'd')
hold on
scatter(GMbestList, m2payloadFractionList_m3p1_m2p4, sz, 'filled', '^')
hold on
scatter(GMbestList, m3antennaLength_m3p1_m2p4, sz, 'filled', 'o')
xlim([10 70])
ylim([0 2.5])
ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'on';
xlabel("GM Best Assumption",'FontSize',18,'FontWeight','bold', 'FontName', 'Times New Roman')
ylabel("USC Best Scoring Plane",'FontSize',18,'FontWeight','bold', 'FontName', 'Times New Roman')
lgd = legend('Span','M2 Payload Fraction','M3 Antenna Length');
lgd.FontSize = 18;
lgd.FontName = 'Times New Roman';

subplot(2,1,2)
hold on
scatter(GMbestList, ARList_m3p1_m2p4, sz, 'filled', 'square')
xlim([10 70])
ylim([10 14])
ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'on';
xlabel("GM Best Assumption",'FontSize',18,'FontWeight','bold', 'FontName', 'Times New Roman')
ylabel("USC Best Scoring Plane",'FontSize',18,'FontWeight','bold', 'FontName', 'Times New Roman')
lgd = legend('AR');
lgd.FontSize = 18;
lgd.FontName = 'Times New Roman';
