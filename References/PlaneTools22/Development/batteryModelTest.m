a = 3.5;
b = -.0334;
c = -.106;
d = .74;
m = 1.4;
n = 2;

nSeries = 1;

Voc = @(s) a + b.*(-log(s)).^m + c.*s + d.*exp(n.*(s-1));
s = .02:.02:1;

Rb = .0012;
Rm = [.014 .02 .05 .08];
capacity = 4.4;

Rt = Rm + Rb;

I = [Voc(s)/Rt(1); Voc(s)/Rt(2); Voc(s)/Rt(3); Voc(s)/Rt(4)];
meanI = mean(transpose(I));
cRate = meanI/capacity;

Vb = nSeries*[Voc(s)-I(1,:)*Rb; Voc(s)-I(2,:)*Rb; Voc(s)-I(3,:)*Rb; Voc(s)-I(4,:)*Rb];

figure(1);
hold on;

plot(1-s, Vb(1,:))
plot(1-s, Vb(2,:))
plot(1-s, Vb(3,:))
plot(1-s, Vb(4,:))

ylim(nSeries*[0 4.5])

title('2.2Ah 1S LiPo battery decay - Beijing Jiaotong University model');
xlabel('State of charge');
ylabel('Battery terminal voltage [V]');
%legend([num2str(meanI(1), 3) ' A'], [num2str(meanI(2), 2) ' A'],...
%   [num2str(meanI(3), 2) ' A'], [num2str(meanI(4), 2) ' A']);
legend([num2str(cRate(1), 2) 'C'], [num2str(cRate(2), 2) 'C'],...
    [num2str(cRate(3), 2) 'C'], [num2str(cRate(4), 2) 'C']);

hold off;