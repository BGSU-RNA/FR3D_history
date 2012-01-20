% zSmallVersion(File) removes some fields from File 

function [File] = zSmallVersion(File)

File.Pair = [];
File.CI   = [];
File.SizeCode = 2;