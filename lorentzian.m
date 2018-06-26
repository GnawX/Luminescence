function [ a ] = lorentzian( x,u,sigma )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

a = sigma^2/((x-u)^2 + sigma^2);
end

