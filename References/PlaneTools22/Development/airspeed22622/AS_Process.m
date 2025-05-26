clear; clc;
close all

%% Flight 1
raw_12 = csvread('Flight 12.csv');

t_ms12 = raw_12(:,1);
t12 = t_ms12/1000;
AS12 = raw_12(:,2);

AS12_P = AS12-8104;
for i=1:length(AS12_P)
    if AS12_P(i)<0
        AS12_P(i)=0;
    end
end
TAS12 = (sqrt(AS12_P)+0.39331)/0.83039;

figure()
plot(t12,TAS12);
title('Flight 1 Airspeed');
xlabel('Time [s]');
ylabel('Airspeed [m/s]');
xlim([390,420]);

csvwrite('Flight 12 TAS.csv',[t12,TAS12]);

%% Flight 2
figure()
plot(t12,TAS12);
title('Flight 2 Airspeed');
xlabel('Time [s]');
ylabel('Airspeed [m/s]');
xlim([1520,1730]);

%% Flight 3
raw_3 = csvread('Flight 3.csv');

t_ms3 = raw_3(:,1);
t3 = t_ms3/1000;
AS3 = raw_3(:,2);

AS3_P = AS3-8104;
for i=1:length(AS3_P)
    if AS3_P(i)<0
        AS3_P(i)=0;
    end
end
TAS3 = (sqrt(AS3_P)+0.39331)/0.83039;

figure()
plot(t3,TAS3);
title('Flight 3 Airspeed');
xlabel('Time [s]');
ylabel('Airspeed [m/s]');
xlim([500,550]);

csvwrite('Flight 3 TAS.csv',[t3,TAS3]);

%% Flight 4
raw_4 = csvread('Flight 4.csv');

t_ms4 = raw_4(:,1);
t4 = t_ms4/1000;
AS4 = raw_4(:,2);

AS4_P = AS4-8104;
for i=1:length(AS4_P)
    if AS4_P(i)<0
        AS4_P(i)=0;
    end
end
TAS4 = (sqrt(AS4_P)+0.39331)/0.83039;

figure()
plot(t4,TAS4);
title('Flight 4 Airspeed');
xlabel('Time [s]');
ylabel('Airspeed [m/s]');
xlim([60,130]);

csvwrite('Flight 4 TAS.csv',[t4,TAS4]);

%% Other Flights
% raw_0 = csvread('Flight 0.csv');
% 
% t_ms0 = raw_0(:,1);
% t0 = t_ms0/1000;
% AS0 = raw_0(:,2);
% 
% figure()
% plot(t0,AS0);
% title('Flight 0 Airspeed');
% xlabel('Time [s]');
% ylabel('Airspeed [ct]');
% xlim([2150,2750]);
