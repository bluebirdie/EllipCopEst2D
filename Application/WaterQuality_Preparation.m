clearvars;

ImportWaterQualityData;

% Temperature
Temp=TOW(:,1);
% Dissolved Nitrogen
DN=TOW(:,end);


%% Check outliers of Temp and DN
figure
subplot(2, 2, 1)
boxplot(Temp)
title("Temperature")
subplot(2, 2, 2)
boxplot(DN)
title("Dissolved nitrogen")
subplot(2, 2, 3)
histogram(Temp, 'Normalization', "pdf")
title("Temperature")
subplot(2, 2, 4)
histogram(DN, 'Normalization', "pdf")
title("Dissolved nitrogen")

%  index of the outliers (only the one differs mostly)
Temp_out=1036;
DN_out=790;

% exclude observations with either outliers in temperature or dissolved
% nitrogen
index_no_out=setdiff(1:length(Temp), [1036 790]);
Temp2=Temp(index_no_out);
DN2=DN(index_no_out);

% boxplot with outliers removed
figure
subplot(2, 2, 1)
boxplot(Temp2)
title("Temperature")
subplot(2, 2, 2)
boxplot(DN2)
title("Dissolved nitrogen")
subplot(2, 2, 3)
histogram(Temp2, 'Normalization', "pdf")
title("Temperature")
subplot(2, 2, 4)
histogram(DN2, 'Normalization', "pdf")
title("Dissolved nitrogen")


%% Check missing values in the data
sum(isnan(Temp2)) 
sum(isnan(DN2))  

% Remove missing observations with either missing temperature or missing
% dissolved nitrogen
index_nomissing=~isnan(Temp2) & ~isnan(DN2);
Temp3=Temp2(index_nomissing);
DN3=DN2(index_nomissing);

% Check missing values in the data
sum(isnan(Temp3)) 
sum(isnan(DN3))  

close all;

% boxplot with outliers and missing values removed
figure
subplot(2, 2, 1)
boxplot(Temp3)
title("Temperature")
subplot(2, 2, 2)
boxplot(DN3)
title("Dissolved nitrogen")
subplot(2, 2, 3)
histogram(Temp3, 'Normalization', "pdf")
title("Temperature")
subplot(2, 2, 4)
histogram(DN3, 'Normalization', "pdf")
title("Dissolved nitrogen")

%% Obtain Pearson correlation 
Temp_DN_data=[Temp3 DN3];
Omega=corr(Temp_DN_data)

%% Obtain Marginal Distributions

pd_temp=fitgmdist(Temp3, 3);
x_values = -10:0.1:50;
y = pdf(pd_temp, x_values');
figure
histogram(Temp3, 'Normalization','pdf', 'BinWidth', 1)
hold on
plot(x_values,y)


pd_DN = fitdist(DN3, 'Kernel');
x_values = 0:0.01:5;
y = pdf(pd_DN,x_values);
figure
histogram(DN3, 'Normalization','pdf')
hold on
plot(x_values,y)

%% Obtain copula data

U_Temp = cdf(pd_temp,Temp3);
figure
histogram(U_Temp, 'Normalization','pdf')

U_DN = cdf(pd_DN,DN3);
figure
histogram(U_DN, 'Normalization','pdf')


figure
histogram2(U_DN, U_Temp, 'Normalization','pdf', 'FaceAlpha', 1, 'FaceColor','#EDB120')


U=[U_Temp U_DN];

save("CopulaData.mat", 'U', 'Omega')


