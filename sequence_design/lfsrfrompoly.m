function c=lfsrfrompoly(inpoly)

% Linear Feedback Shift Register for generating pseudorandom binary sequences
%
% Written by Ben Kravitz (ben.kravitz@pnnl.gov or ben.kravitz.work@gmail.com)
% 11 May 2016
%
% Inputs:
%   inpoly = the polynomial (excluding the 0th term) used to generate the sequence
%      example:  [12 10 2 1] = x^12 + x^10 + x^2 + x + 1
%
% Will generate a sequence of 2^s-1 (maximum length possible)
% Initial state is all ones
% Uses the fibonacci method
% You are responsible for ensuring the polynomial you enter is primitive.  The
%   script will not do this for you.

s=max(inpoly);

c=ones(1,s);
while length(c)<2^s-1;
    worksq=c(1:s);
    j=worksq(inpoly(1));
    for k=2:length(inpoly);
        j=xor(j,worksq(inpoly(k)));
    end
    c=[j c];
end
