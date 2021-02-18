%% for generate the interface defination file 
%% of ClassCalculateA20InOneElement.hpp 
%% and ClassCalculateA20InOneElement.cpp
%% labrary

clear
clc
productPath='./';
hppFileMy='ClassCalculateA20InOneElement.hpp';
cppFileMy='ClassCalculateA20InOneElement.cpp';

clibgen.generateLibraryDefinition(fullfile(productPath,hppFileMy),...
"SupportingSourceFiles",fullfile(productPath,cppFileMy),...
"IncludePath",productPath,...
"ReturnCArrays",false,...
'Verbose',true);
