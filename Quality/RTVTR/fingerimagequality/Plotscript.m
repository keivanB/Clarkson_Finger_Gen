%plotting sample script for paper

clear all
close all
clc
Feature ='ridgevalleythiknessratio';
File = '/home/citer/Desktop/Fingerprint_Age/Clarkson_Age_Estimation/Processed Images/Dataset/Resized/ALLinone/034_120033_Image_00.bmp';
X = 6;
Y = 12;
[Q_RVU, Q_Ratios, Q_RVU_local, RTVTR_Ratios, Q_profile, Q_ridval, Q_change, Q_Blocks, Q_R_Blocks,Q_WLC, Q_RC] = RTVTR_Interface(Feature,File);
figure; imshow(imresize(Q_Blocks{X, Y},20));
%figure; imshow(Q_R_Blocks{4, 10});
figure; plot(Q_profile{X, Y});
%rirges 1 , valleys 0
figure; bar(Q_ridval{X, Y});
figure; bar(Q_change{X, Y});
RVU_vector = reshape(Q_RVU_local,1,[]);
RVU_vector = sort(RVU_vector(~isnan(RVU_vector)));
RVU_vector = RVU_vector(RVU_vector>0);