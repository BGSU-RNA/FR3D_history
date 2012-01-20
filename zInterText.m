% zInterText returns a text string telling the interaction between NT1 and NT2

function [Text] = zInterText(NT1,NT2,Class,display)   
    %display = 1 means short annotation, any other value means long annotation

Num1 = '    ';
Num1((5-length(NT1.Number)):4) = NT1.Number;  %make Num1, Num2 have same length

Num2 = '    ';
Num2((5-length(NT2.Number)):4) = NT2.Number;

Base1 = [NT1.Base ' ' Num1 ' (' NT1.Chain ') ']; % text for Base 1
Base2 = [NT2.Base ' ' Num2 ' (' NT2.Chain ') ']; % text for Base 2

if NT1.Base < NT2.Base                           
  Low  = Base1;
  High = Base2; 
elseif NT2.Base < NT1.Base,
  Low  = Base2;
  High = Base1;
elseif Num1 < Num2,                 % if NT1.Base = NT2.Base, look again
  Low  = Base1;
  High = Base2; 
elseif Num1 < Num2,
  Low  = Base2;
  High = Base1;
elseif NT1.Chain <= NT2.Chain,      % if Bases and numbers are the same
  Low  = Base1;
  High = Base2; 
else
  Low  = Base2;
  High = Base1;
end

if Class > 0,
  Dom = Low;                           % base using the dominant edge
  Sec = High;                          % base using the dominant edge
else
  Dom = High;                          % base using the dominant edge
  Sec = Low;                           % base using the secondary edge
end

% List base using dominant edge first
if display == 1,
    switch abs(Class)
        case   1,    Text = ['cWW'];
        case   2,    Text = ['tWW'];
        case   3,    Text = ['cWH'];
        case   4,    Text = ['tWH'];
        case   5,    Text = ['cWS'];
        case   6,    Text = ['tWS'];
        case   7,    Text = ['cHH'];
        case   8,    Text = ['tHH'];
        case   9,    Text = ['cHS'];
        case  10,    Text = ['tHS'];
        case  11,    Text = ['cSS'];
        case  12,    Text = ['tSS'];
        case  13,    Text = ['bif'];
        case  25,    Text = ['part of motif'];
        otherwise,   Text = ['-'];  % add cases later
    end
else 
    switch abs(Class)
        case   1,    Text = ['Cis   ' Dom  'WC - ' Sec 'WC'];
        case   2,    Text = ['Trans ' Dom  'WC - ' Sec 'WC'];
        case   3,    Text = ['Cis   ' Dom  'WC - ' Sec 'H '];
        case   4,    Text = ['Trans ' Dom  'WC - ' Sec 'H '];
        case   5,    Text = ['Cis   ' Dom  'WC - ' Sec 'S '];
        case   6,    Text = ['Trans ' Dom  'WC - ' Sec 'S '];
        case   7,    Text = ['Cis   ' Dom  'H  - ' Sec 'H '];
        case   8,    Text = ['Trans ' Dom  'H  - ' Sec 'H '];
        case   9,    Text = ['Cis   ' Dom  'H  - ' Sec 'S '];
        case  10,    Text = ['Trans ' Dom  'H  - ' Sec 'S '];
        case  11,    Text = ['Cis   ' Dom  'S  - ' Sec 'S '];
        case  12,    Text = ['Trans ' Dom  'S  - ' Sec 'S '];
        otherwise,   Text = ['      ' Dom  '     ' Sec '  '];  % add cases later
    end
end