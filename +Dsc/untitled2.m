clear all
clc
close all
%%
x = [-10:0.1:10];
y =  (1 + exp(-  x )).^-1;
plot(x,y,'linewidth',2)
title('y =  (1 + exp(-  x )).^-1','fontsize',18)
set(gca,'fontsize',18)
figure
x1= [0:0.1:20];
y1 = (1 + exp(-  (x1-10) )).^-1;
plot(x1,y1,'linewidth',2)
title('y = (1 + exp(-  (x-10) )).^-1','fontsize',18)
set(gca,'fontsize',18)
hold on
sensitivity = 2;
bias = 5;
y2 = (1 + exp(-  (bias+sensitivity*(x1-10)) )).^-1;
plot(x1,y2,'linewidth',2)
%%
% rewrite the generative model
N = 1000;
x =  20*rand(N,1);

sensitivity = 2;
bias = 5;
x_mean = 10;
p = (1 + exp(-  (bias+sensitivity*(x-10)) )).^-1;

choice = zeros(N,1);
rnds = rand(N,1);
choice(rnds < p ) = 1;
figure
plot(x,p,'.')
hold on
plot(x,choice,'*')
% run this multiple times
% see if empirical p matches theorical p
nRun = 200;
choiceMulti = zeros(N,nRun);

for i = 1:nRun
    rnds = rand(N,1);
    choiceMulti(rnds < p,i) = 1;
end
figure
plot(x,p,'.')
hold on
plot(x,mean(choiceMulti,2),'*');