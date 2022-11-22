clc
clear
close all
%% Load Data %%
load('lefthand.mat');
load('righthand.mat');
load('left.mat');

leftHand = table2array(lefthand);
rightHand = table2array(righthand);
left = table2array(left);

sampleFre = 20;
deltT = 1/sampleFre;
tl = (1:length(leftHand)) * deltT;
tr = (1:length(rightHand)) * deltT;
tle = (1:length(left)) * deltT;
%% Data Segmentation %%
% left hand
axl = leftHand(:,2);
% segmentation
axl_1 = axl(1:182);
axl_2 = axl(182:412);
axl_3 = axl(412:662);
axl_4 = axl(662:922);
axl_5 = axl(922:1138);
axl = {axl_1,axl_2,axl_3,axl_4,axl_5};

% right hand
axr = rightHand(:,2);
% segmentation
axr_1 = axr(1:148);
axr_2 = axr(148:310);
axr_3 = axr(310:495);
axr_4 = axr(495:674);
axr_5 = axr(674:857);
axr = {axr_1,axr_2,axr_3,axr_4,axr_5};

% left hand 20 
al = left(1:4000,2);
al1 = al(1:203);
al2 = al(225:398);
al3 = al(430:589);
al4 = al(633:753);
al5 = al(829:959);
al6 = al(1005:1137);
al7 = al(1195:1321);
al8 = al(1388:1532);
al9 = al(1580:1710);
al10 = al(1772:1892);
al11 = al(1984:2119);
al12 = al(2171:2321);
al13 = al(2354:2509);
al14 = al(2533:2668);
al15 = al(2709:2855);
al16 = al(2895:3027);
al17 = al(3063:3194);
al18 = al(3219:3319);
al19 = al(3372:3486);
al = {al1,al2,al3,al4,al5,al6,al7,al8,al9,al10,al11,al12,al13,al14,al15,al16,al17,al18,al19};

%% Window Length Statistics %%
leftwindow = [length(axl_1),length(axl_2),length(axl_3),length(axl_4),length(axl_5)];
rightwindow = [length(axr_1),length(axr_2),length(axr_3),length(axr_4),length(axr_5)];
leftwindow_ave = mean(leftwindow);  leftwindow_std = std(leftwindow);
rightwindo_ave = mean(rightwindow);  rightwindow_std = std(rightwindow);
%% Filter Application %%
for i = 1:length(axl)
    axl_f{i} = preprocess_sgolay(axl{i});
    axr_f{i} = preprocess_sgolay(axr{i});
end
%% Z Normalisation %%
for i = 1:length(axl_f)
    axl_n{i} = zscore(axl_f{i});
    axr_n{i} = zscore(axr_f{i});
end
%% normalisation and filter
for i = 1:length(al)
    aln{i} = preprocess_sgolay_n(al{i});
end

for i = 1:length(al)
    alf{i} = preprocess_sgolay(al{i});
end
%% DTW %%
% left hand 
[dl1,pl1,ql1] = dtw(axl_n{1}, axl_n{2});
[dl2,pl2,ql2] = dtw(axl_n{1}, axl_n{3});
[dl3,pl3,ql3] = dtw(axl_n{1}, axl_n{4});
[dl4,pl4,ql4] = dtw(axl_n{1}, axl_n{5});

[dl5,pl5,ql5] = dtw(axl_n{2}, axl_n{3});
[dl6,pl6,ql6] = dtw(axl_n{2}, axl_n{4});
[dl7,pl7,ql7] = dtw(axl_n{2}, axl_n{5});

[dl8,pl8,ql8] = dtw(axl_n{3}, axl_n{4});
[dl9,pl9,ql9] = dtw(axl_n{3}, axl_n{5});

[dl10,pl10,ql10] = dtw(axl_n{4}, axl_n{5});

dl = [dl1,dl2,dl3,dl4,dl5,dl6,dl7,dl8,dl9,dl10];
dl_ave = mean(dl);   dl_std = std(dl);

% right hand
[dr1,pr1,qr1] = dtw(axr_n{1}, axr_n{2});
[dr2,pr2,qr2] = dtw(axr_n{1}, axr_n{3});
[dr3,pr3,qr3] = dtw(axr_n{1}, axr_n{4});
[dr4,pr4,qr4] = dtw(axr_n{1}, axr_n{5});

[dr5,pr5,qr5] = dtw(axr_n{2}, axr_n{3});
[dr6,pr6,qr6] = dtw(axr_n{2}, axr_n{4});
[dr7,pr7,qr7] = dtw(axr_n{2}, axr_n{5});

[dr8,pr8,qr8] = dtw(axr_n{3}, axr_n{4});
[dr9,pr9,qr9] = dtw(axr_n{3}, axr_n{5});

[dr10,pr10,qr10] = dtw(axr_n{4}, axr_n{5});

dr = [dr1,dr2,dr3,dr4,dr5,dr6,dr7,dr8,dr9,dr10];
dr_ave = mean(dr);   dr_std = std(dr);

% left-right hand
for i = 1:length(aln)
    [d(i),p{i},q{i}] = dtw(axr_n{2},aln{i});
end
%% Sample Entropy
SE_temp = SampEntropy(2, 0.2*std(axr_n{2}), axr_n{2});

for i = 1:length(alf)
    SE(i) = SampEntropy(2, 0.2*std(alf{i}),alf{i});
end
%% Spectral Entropy %%
for i = 1:length(alf)
    E(i) = SpecEntropy(alf{i});
end
%% Visualisation
% DTW
x = 1:length(d);
[p,S] = polyfit(x,d,1);
[Rd,Pd] = corrcoef(x,d);
[D1,delta] = polyval(p,x,S);
figure
plot(x,d,'bo');
hold on 
plot(x,D1,'r-');
plot(x,D1+2*delta,'m--',x,D1-2*delta,'m--')
str1 = {'Correlation coefficient: -0.63','p=0.0038'};
dim1 = [0.15 0.15 0.1 0.1];
annotation('textbox',dim1,'String',str1,'FitBoxToText','on');
legend('Data','Linear Fit','95% Prediction Interval')
title('Linear Fit of Data with 95% Prediction Interval-DTW')

% Sample Entropy
[ps,Ss] = polyfit(x,SE,1);
[Rs,Ps] = corrcoef(x,SE);
[SE1,deltas] = polyval(ps,x,Ss);
figure
plot(x,SE,'bo');
hold on 
plot(x,SE1,'r-')
plot(x,SE1+2*deltas,'m--',x,SE1-2*deltas,'m--')
str2 = {'Correlation coefficient: -0.5582','p=0.0130'};
annotation('textbox',dim1,'String',str2,'FitBoxToText','on');
legend('Data','Linear Fit','95% Prediction Interval')
title('Linear Fit of Data with 95% Prediction Interval-Sample Entropy')

% DTW vs Sample Entropy
figure
plot(d,SE,'bo');
title('DTW against Sample Entropy')
xlabel('DTW')
ylabel('Sample Entropy')

[Rds,Pds] = corrcoef(d,SE);



