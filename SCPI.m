clc,clear,close all
load 'data.mat'

t=out1.simout_debit.Time(81:196)-8;
y1=out1.simout_debit.Data(81:196)-1;
y2=out2.simout_debit.Data(81:196)-1;
y3=out3.simout_debit.Data(81:196)-1;

plot(t,y1)
hold on
plot(t,y2)
plot(t,y3)
xlabel("Timp (s)")
legend("y_1","y_2","y_3",'Location','southeast')
title("Raspunsurile indiciale")

%% filtru butter
[b,a]=butter(2,0.3,'low');
y1_f=filtfilt(b,a,y1);
y2_f=filtfilt(b,a,y2);
y3_f=filtfilt(b,a,y3);

figure
plot(t,y1_f)
hold on
plot(t,y2_f)
plot(t,y3_f)
xlabel("Timp (s)")
legend("y_{1f}","y_{2f}","y_{3f}",'Location','southeast')
title("Raspunsurile indiciale filtrate")

%% mediere + functie de transfer ord1
y_f=(y1_f+y2_f+y3_f)/3;
figure
plot(t,y_f)
hold on

info=stepinfo(y_f,t,'SettlingTimeThreshold',0.05)
K=info.SettlingMax/160
T=info.TransientTime/3
s=tf('s');
H=tf(K,[T 1])*exp(-0.2*s)

u=160.*double(t>=0);
y_sim=lsim(H,u,t);
plot(t,y_sim)
xlabel("Timp (s)")
legend("y_f","y_{sim}",'Location','southeast')
title("Raspunsul indicial masurat, filtrat si mediat si cel simulat")

%% regulator P1
Kr=T/K;
Ti=T;
Td=0;

load data_regulator.mat
figure
plot(out.simout_debit)
title("Validare regulator P1")

% robustete
H_pid=Kr*(1+1/(Ti*s));
L=H_pid*H/(1+H_pid*(tf(K,[T 1])-H));
S=feedback(1,L)
figure
bode(S)
deltaM=1/norm(pade(S,2),'inf')

%% regulator P2
H2=0.06/(0.206*s+1)/(0.24*s+1)/(12*s+1)*exp(-0.2*s)

Kr2=20.4;
Ti2=12.24;
Td2=0.235;

%% cascada

% GT interior 
zeta=0.95;
tt=2;
wn=4/zeta/tt;
GT=tf(wn^2,[1 2*wn*zeta wn^2])

load date_cascada.mat
figure
plot(out.simout_nivel)
hold on
plot(out.comanda_nivel)
ylabel("")
legend("y","u",'Location','southeast')
title("Validare cascada-nivel")

figure
plot(out.simout_debit)
hold on
plot(out.comanda_debit)
ylabel("")
legend("y","u",'Location','southeast')
title("Validare cascada-debit")