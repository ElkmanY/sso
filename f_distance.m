function [ Dist ] = f_distance( subject, object  )    
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%global N D
Dist=sqrt(sum((subject-object).^2,2));











end

