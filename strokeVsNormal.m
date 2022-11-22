clc
clear
close all

%% Duke data %%
sampleRate = 20; % Hz
deltT = 1 / sampleRate;
dataNormal = csvread('normal.csv', 1, 0, [1, 0, 551, 3]);
timeNormal = dataNormal(:,1) * deltT;
dataStroke = csvread('stroke.csv', 1, 0, [1, 0, 1550, 3]);
timeStroke = dataStroke(:,1) * deltT;

zNormal = dataNormal(:,4);
xNormal = dataNormal(:,2);
yNormal = dataNormal(:,3);

zStroke = dataStroke(:,4);
xStroke = dataStroke(:,2);
yStroke = dataStroke(:,3);

%% Rayna data %%
dataRN = csvread('RaynaNormal.csv', 1, 0, [1, 0, 729, 3]);
timeRN = dataRN(:,1) * deltT;
dataStrokeR = csvread('strokeRayna.csv',1,0,[1,0,922,3]);
timeStrokeR = dataStrokeR(:,1) * deltT;

zRNormal = dataRN(:,4);
zStrokeR = dataStrokeR(:,4);

%% Kim data %%
dataKimN = csvread('NormalKim.csv',1,0,[1,0,561,3]);
dataKimStroke = csvread('StrokeKim.csv',1,0,[1,0,763,3]);
timeKimN = dataKimN(:,1) * deltT;
timeKimStroke = dataKimStroke(:,1) * deltT;

zNormalKim = dataKimN(:,4);
zStrokeKim = dataKimStroke(:,4);

%% Raw data plot %%
% figure
% subplot(2,1,1)
% plot(timeNormal, zNormal);
% title('Normal Reaching Shoulder Movement Acceleration (Duke)')
% xlabel('Time (s)');
% ylabel('Z-axis Acceleration (g)');
% grid on
% 
% subplot(2,1,2)
% plot(timeStroke, zStroke);
% title('Stroke Reaching Shoulder Movement Acceleration')
% xlabel('Time (s)');
% ylabel('Z-axis Acceleration (g)');
% grid on
% 
% figure
% subplot(2,1,1)
% plot(timeRN, zRNormal);
% title('Normal Reaching Shoulder Movement Acceleration (Rayna)')
% xlabel('Time (s)');
% ylabel('Z-axis Acceleration (g)');
% grid on
% 
% subplot(2,1,2)
% plot(timeStrokeR, zStrokeR);
% title('Stroke Reaching Shoulder Movement Acceleration (Rayna)')
% xlabel('Time (s)');
% ylabel('Z-axis Acceleration (g)');
% grid on
% 
% figure
% subplot(2,1,1)
% plot(timeKimN, zNormalKim);
% title('Normal Reaching Shoulder Movement Acceleration (Kim)')
% xlabel('Time (s)');
% ylabel('Z-axis Acceleration (g)');
% grid on
% 
% subplot(2,1,2)
% plot(timeKimStroke,zStrokeKim)
% title('Stroke Reaching Shoulder Movement Acceleration (Kim)')
% xlabel('Time (s)');
% ylabel('Z-axis Acceleration (g)');
% grid on

%% Data filtering %%
zNormalF = preprocess(zNormal);
zRNormalF = preprocess(zRNormal);
zStrokeF = preprocess(zStroke);
zStrokeRF = preprocess(zStrokeR);
zNormalKimF = preprocess(zNormalKim);
zStrokeKimF = preprocess(zStrokeKim);

figure
subplot(2,3,1)
plot(timeNormal, zNormalF);
title('Normal Acceleration (Duke)')
xlabel('Time (s)');
ylabel('Z-axis Acceleration (g)');
grid on
subplot(2,3,4)
plot(timeStroke, zStrokeF);
title('Stoke Acceleration (Duke)')
xlabel('Time (s)');
ylabel('Z-axis Acceleration (g)');
grid on

subplot(2,3,2)
plot(timeRN, zRNormalF);
title('Normal Acceleration (Rayna)')
xlabel('Time (s)');
ylabel('Z-axis Acceleration (g)');
grid on
subplot(2,3,5)
plot(timeStrokeR, zStrokeRF);
title('Stroke Acceleration (Rayna)')
xlabel('Time (s)');
ylabel('Z-axis Acceleration (g)');
grid on

subplot(2,3,3)
plot(timeKimN, zNormalKimF);
title('Normal Acceleration (Kim)')
xlabel('Time (s)');
ylabel('Z-axis Acceleration (g)');
grid on
subplot(2,3,6)
plot(timeKimStroke, zStrokeKimF);
title('Stroke Acceleration (Kim)')
xlabel('Time (s)');
ylabel('Z-axis Acceleration (g)');
grid on

sgtitle('Filtered Acceleration (Gaussian Filter with 10 Samples Window) From Three Subjects (1 Male and 2 Females)')

%% Data normalisation %% 
zNormalN = zscore(zNormalF);
zRNormalN = zscore(zRNormalF);
zStrokeN = zscore(zStrokeF);
zStrokeRN = zscore(zStrokeRF);
zStrokeKimN = zscore(zStrokeKimF);
zNormalKimN = zscore(zNormalKimF);

csvwrite('zNormalDuke.csv',zNormalN);
csvwrite('zNormalRay.csv',zRNormalN);
csvwrite('zNormalKim.csv',zNormalKimN);
csvwrite('zStrokeDuke.csv',zStrokeN);
csvwrite('zStrokeRay.csv',zStrokeRN);
csvwrite('zStrokeKim.csv',zStrokeKimN);

%% DTW calculation %% 
figure 
dtw(zNormalN, zRNormalN);
xlabel('Warped Dukes healthy data and Raynas healthy data');
[d1,i11,i12] = dtw(zNormalN, zRNormalN);
r1 = max(round(0.25*max(length(zNormalN),length(zRNormalN))), abs(length(zNormalN)-length(zRNormalN)));
d1Regularised = dtw(zNormalN, zRNormalN, r1);

figure
dtw(zNormalN, zNormalKimN);
xlabel('Warped Dukes healthy data and Kims healthy data');
[d2,i21,i22] = dtw(zNormalN, zNormalKimN);
r2 = max(round(0.25*max(length(zNormalN),length(zNormalKimN))), abs(length(zNormalN)-length(zNormalKimN)));
d2Regularised = dtw(zNormalN, zNormalKimN, r2);

figure
dtw(zNormalKimN, zRNormalN);
xlabel('Warped Raynas healthy data and Kims healthy data');
[d3,i31,i32] = dtw(zNormalKimN, zRNormalN);
r3 = max(round(0.25*max(length(zNormalKimN),length(zRNormalN))), abs(length(zNormalKimN)-length(zRNormalN)));
d3Regularised = dtw(zNormalKimN, zRNormalN, r3);

figure
dtw(zNormalN, zStrokeN);
xlabel('Warped Dukes healthy data and Dukes stroke-like data');
[d4,i41,i42] = dtw(zNormalN, zStrokeN);
r4 = max(round(0.25*max(length(zNormalN),length(zStrokeN))), abs(length(zNormalN)-length(zStrokeN)));
d4Regularised = dtw(zNormalN, zStrokeN, r4);

figure
dtw(zRNormalN, zStrokeRN);
xlabel('Warped Raynas healthy data and Raynas stroke-like data');
[d5,i51,i52] = dtw(zRNormalN, zStrokeRN);
r5 = max(round(0.25*max(length(zRNormalN),length(zStrokeRN))), abs(length(zStrokeRN)-length(zRNormalN)));
d5Regularised = dtw(zStrokeRN, zRNormalN, r5);

figure
dtw(zNormalKimN, zStrokeKimN);
xlabel('Warped Kims healthy data and Kims stroke-like data');
[d6,i61,i62] = dtw(zNormalKimN, zStrokeKimN);
r6 = max(round(0.25*max(length(zNormalKimN),length(zStrokeKimN))), abs(length(zNormalKimN)-length(zStrokeKimN)));
d6Regularised = dtw(zNormalKimN, zStrokeKimN, r6);

dHealthy = [d1,d2,d3];
dStroke = [d4,d5,d6];
D = [dHealthy',dStroke'];

figure
boxplot(D, 'Labels',{'Healthy','Stroke'});
dim1 = [0.2 0.5 0.3 0.3];
dim2 = [0.5 0.1 0.3 0.3];
str1 = {'Mean value: 73.0325','Standard deviation: 29.7221'};
str2 = {'Mean value: 147.5451','Standard deviation: 26.2721'};
annotation('textbox',dim1,'String',str1,'FitBoxToText','on');
annotation('textbox',dim2,'String',str2,'FitBoxToText','on');
title('Healthy and Stroke DST Comparison');
ylabel('Minimum distance');

dHealthy_ave = sum(dHealthy)/3;
dHealthy_std = std(dHealthy);
dStroke_ave = sum(dStroke)/3;
dStroke_std = std(dStroke);
fprintf('The average dynamic distance of healthy group is %f with a standard deviation of %f \n',dHealthy_ave, dHealthy_std);
fprintf('The average dynamic distance of stroke group is %f with a standard deviation of %f \n',dStroke_ave, dStroke_std);

%% visualisation %% 
for i = 1:length(zNormalN)
    for j = 1:length(zRNormalN)
        dis1(i,j) = abs(zNormalN(i)-zRNormalN(j));
    end
end
figure
subplot(2,3,1)
p1 = pcolor(dis1);
set(p1, 'EdgeColor','none');
colorbar
hold on
plot(i12, i11, 'r-','linewidth',3);
title('Duke vs Rayna');

for i = 1:length(zNormalN)
    for j = 1:length(zNormalKimN)
        dis2(i,j) = abs(zNormalN(i) - zNormalKimN(j));
    end
end
subplot(2,3,2)
p2 = pcolor(dis2);
set(p2, 'EdgeColor','none');
colorbar
hold on
plot(i22, i21, 'r-','linewidth',3);
title('Duke vs Kim');

for i = 1:length(zRNormalN)
    for j = 1:length(zNormalKimN)
        dis3(i,j) = abs(zRNormalN(i) - zNormalKimN(j));
    end
end
subplot(2,3,3)
p3 = pcolor(dis3);
set(p3, 'EdgeColor','none');
colorbar
hold on
plot(i31, i32, 'r-','linewidth',3);
title('Rayna vs Kim');

for i = 1:length(zNormalN)
    for j = 1:length(zStrokeN)
        dis4(i,j) = abs(zNormalN(i) - zStrokeN(j));
    end
end
subplot(2,3,4)
p4 = pcolor(dis4);
set(p4, 'EdgeColor','none');
colorbar
hold on
plot(i42, i41, 'r-','linewidth',3);
title('Duke healthy vs Duke stroke');

for i = 1:length(zRNormalN)
    for j = 1:length(zStrokeRN)
        dis5(i,j) = abs(zRNormalN(i) - zStrokeRN(j));
    end
end
subplot(2,3,5)
p5 = pcolor(dis5);
set(p5, 'EdgeColor','none');
colorbar
hold on
plot(i52, i51, 'r-','linewidth',3);
title('Rayna healthy vs Rayna stroke');

for i = 1:length(zNormalKimN)
    for j = 1:length(zStrokeKimN)
        dis6(i,j) = abs(zNormalKimN(i) - zStrokeKimN(j));
    end
end
subplot(2,3,6)
p6 = pcolor(dis6);
set(p6, 'EdgeColor','none');
colorbar
hold on
plot(i62, i61, 'r-','linewidth',3);
title('Kim healthy vs Kim stroke');
sgtitle('Minimum Distance Among the Distance Matrix');

%% Cross-correlation %%
% [c1,lags1]=xcorr(zNormalF,zRNormalF);
% c1_max = max(c1);
% lags1_c1_max = lags1(find(c1==c1_max));
% 
% [c2,lags2]=xcorr(zNormalF,zNormalKimF);
% c2_max = max(c2);
% lags2_c2_max = lags2(find(c2==c2_max));
% 
% [c3,lags3]=xcorr(zNormalKimF,zRNormalF);
% c3_max = max(c3);
% lags3_c3_max = lags3(find(c3==c3_max));
% 
% [c4,lags4]=xcorr(zNormalF,zStrokeF);
% c4_max = max(c4);
% lags4_c4_max = lags4(find(c4==c4_max));
% 
% [c5,lags5]=xcorr(zNormalKimF,zStrokeKimF);
% c5_max = max(c5);
% lags5_c5_max = lags5(find(c5==c5_max));
% 
% [c6,lags6]=xcorr(zRNormalF,zStrokeRF);
% c6_max = max(c6);
% lags6_c6_max = lags6(find(c6==c6_max));
% 
% s1 = dukescore(c1_max, lags1_c1_max);
% s2 = dukescore(c2_max, lags2_c2_max);
% s3 = dukescore(c3_max, lags3_c3_max);
% s4 = dukescore(c4_max, lags4_c4_max);
% s5 = dukescore(c5_max, lags5_c5_max);
% s6 = dukescore(c6_max, lags6_c6_max);

%% MJC score from python %%
a1 = 0.026873229400000007;
a2 = 0.004414681400000003;
a3 = 0.002328442599999996;

a4 = 3.485e-6;
a5 = 0.004916746399999986;
a6 = 0.007368880800000003;

aHealthy = [a1,a2,a3];
aStroke = [a4,a5,a6];
A = [aHealthy',aStroke'];

aHealthy_ave = mean(aHealthy);
aHealthy_std = std(aHealthy);
aStroke_ave = mean(aStroke);
aStroke_std = std(aStroke);

figure
subplot(1,2,1)
plot(aHealthy, dHealthy,'*-');
xlabel('MJC dissimilarity');
ylabel('DTW minimum distance');
title('Correlation between MJC and DTW (Healthy data)');
subplot(1,2,2)
plot(aStroke, dStroke,'+-');
xlabel('MJC dissimilarity');
ylabel('DTW minimum distance');
title('Correlation between MJC and DTW (Stroke data)');

figure
boxplot(D, 'Labels',{'Healthy','Stroke'});
dim1 = [0.2 0.5 0.3 0.3];
dim2 = [0.5 0.1 0.3 0.3];
str1 = {'Mean value: 0.0112','Standard deviation: 0.0136'};
str2 = {'Mean value: 0.0041','Standard deviation: 0.0038'};
annotation('textbox',dim1,'String',str1,'FitBoxToText','on');
annotation('textbox',dim2,'String',str2,'FitBoxToText','on');
title('Healthy and Stroke MJC Comparison');
ylabel('Dissimilarity');

%% Spectral Entropy %%
xt1 = timetable(seconds(timeNormal),zNormalF);
e1 = pentropy(xt1).SE / log2(length(zNormalF));
E1 = sum(e1);

xt2 = timetable(seconds(timeStroke),zStrokeF);
e2 = pentropy(xt2).SE / log2(length(zStrokeF));
E2 = sum(e2);

xt3 = timetable(seconds(timeRN),zRNormalF);
e3 = pentropy(xt3).SE / log2(length(zRNormalF));
E3 = sum(e3);

xt4 = timetable(seconds(timeStrokeR),zStrokeRF);
e4 = pentropy(xt4).SE / log2(length(zStrokeRF));
E4 = sum(e4);

xt5 = timetable(seconds(timeKimN),zNormalKimF);
e5 = pentropy(xt5).SE / log2(length(zNormalKimF));
E5 = sum(e5);

xt6 = timetable(seconds(timeKimStroke),zStrokeKimF);
e6 = pentropy(xt6).SE / log2(length(zStrokeKimF));
E6 = sum(e6);




