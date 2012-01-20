% zInteractionName(Class) displays the table name for the given Class

function [Name] = zInteractionName(Class)

switch Class,
  case  1, Name = 'Cis Watson-Crick/Watson-Crick';
  case  2, Name = 'Trans Watson-Crick/Watson-Crick';
  case  3, Name = 'Cis Watson-Crick/Hoogsteen';
  case  4, Name = 'Trans Watson-Crick/Hoogsteen';
  case  5, Name = 'Cis Watson-Crick/Sugar Edge';
  case  6, Name = 'Trans Watson-Crick/Sugar Edge';
  case  7, Name = 'Cis Hoogsteen/Hoogsteen';
  case  8, Name = 'Trans Hoogsteen/Hoogsteen';
  case  9, Name = 'Cis Hoogsteen/Sugar Edge';
  case 10, Name = 'Trans Hoogsteen/Sugar Edge';
  case 11, Name = 'Cis Sugar Edge/Sugar Edge';
  case 12, Name = 'Trans Sugar Edge/Sugar Edge';
end
