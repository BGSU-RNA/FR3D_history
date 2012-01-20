
clear Message

Message{1} = 'Set limits on the difference between nucleotide numbers using the boxes below the diagonal. (Actually, what is used is the difference between the index of nucleotides in the file, not NDB nucleotide number.) To put an upper limit on the difference, type something like <5 or <=5.  To put a lower limit on the difference, type something like >5 or >=5.  To put both limits at once, type something like >5 <=12.  To insist that the nucleotide in the given row have a lower nucleotide number than the nucleotide in the given column, type <, separated by a space from other specifications.  For greater, type >.';

msgbox(Message,'Distance help');

break

fprintf('Set limits on the difference between nucleotide numbers.\n');
fprintf('(Actually, what is used is the difference between the index of\n');
fprintf('nucleotides in the file, not NDB nucleotide number.)\n');
fprintf('To put an upper limit on the difference, type something like <5 or <=5\n');
fprintf('To put a lower limit on the difference, type something like >5 or >=5\n');
fprintf('To put both limits at once, type something like >5 <=12\n');
