clc, clear all; close all;
ds = datastore(fullfile('C:\Users\user\Desktop\Data\atrw_pose_train\train'),...
'IncludeSubfolders', true,'Type', 'image');
fname = 'C:\Users\user\Desktop\Data\atrw_anno_pose_train\keypoint_train.json'; 
fid = fopen(fname); 
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid); 
val = jsondecode(str);
