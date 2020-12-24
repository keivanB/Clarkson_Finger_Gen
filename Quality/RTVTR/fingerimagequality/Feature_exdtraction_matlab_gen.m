%plotting sample script for paper

clear all
close all
clc

d = '/home/citer/Desktop/GANS/Generated_images_512_BMP/Gen_512_BMP/gray/';
fet = '/home/citer/Desktop/GANS/Generated_images_512_BMP/RTVTR/Generated_RTVTR/';
f = dir(d);
files = size(f);
fn = files(1);
Feature ='ridgevalleythiknessratio';

for i=3:fn
    NAME = f(i).name;
    File = append(d,f(i).name);
    [Q_RVU, Q_Ratios, Q_RVU_local, RTVTR_Ratios, Q_profile, Q_ridval, Q_change, Q_Blocks, Q_R_Blocks,Q_WLC, Q_RC] = RTVTR_Interface(Feature, File);
    %[Q_RVU, Q_Ratios, Q_RVU_local, ~, ~, ~, ~, ~, ~,Q_WLC, Q_RC] = RTVTR_Interface(Feature, File);
    %clearvars Q_RVU Q_Blocks Q_change RTVTR_Ratios Q_profile Q_ridval Q_change Q_Blocks Q_R_Blocks
    save(append(fet,NAME,".mat"), 'Q_RVU', 'Q_RVU_local', 'Q_Ratios', 'Q_WLC', 'Q_RC', 'Q_profile')
end

