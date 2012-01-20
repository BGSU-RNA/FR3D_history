
function [Sub] = zSubFile(File,i)

Sub = File;

Sub.NT = File.NT(i);
Sub.Edge = File.Edge(i,i);
Sub.BasePhosphate = File.BasePhosphate(i,i);
Sub.Range = File.Range(i,i);
Sub.NumNT = length(i);
