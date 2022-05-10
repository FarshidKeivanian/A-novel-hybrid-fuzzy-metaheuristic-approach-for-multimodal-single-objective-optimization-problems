clc;
clear;
close all;

load('BestResults.mat');
load('MeanResults.mat');
load('MedianResults.mat');
load('SDResults.mat');

disp(['Best = ',num2str(BestResults)]);
disp(['Mean = ',num2str(MeanResults)]);
disp(['Median = ',num2str(MedianResults)]);
disp(['SD = ',num2str(SDResults)]);