function [y_new,Y_NEW] = matchphases(x,y)
% MATCHPHASES alters the phase of the output, y, to match the phase of the 
% input, x.
%
% written by: Bethany Sutherland
% last edited: 1/22/2019
%
% Inputs:
%   x       = input signal
%   y       = output signal (phase of y will be set to match x)
%
% Outputs:
%
%   y_new   = magnitude squared coherence between input and output signals
%   spec    = the output in the frequency space


szy = size(y);

% take fourier transforms

if size(szy,2) == 3 % make sure you are doing transform along time dimension
    L = min([szy(3), length(x)]);
    Y = fft(y,L,3);
    X = fft(x',L);
else
    L = length(y);
    Y = fft(y,L);
    X = fft(x',L); 
end

% determine magnitudes
Xmag = abs(X);
Ymag = abs(Y);

% determine angles
PhiX = angle(X);

y_new = NaN(size(y));

if size(szy,2) == 3 % make sure you are doing transform along time dimension
    y_new = y_new(:,:,1:L);
else
    y_new = y_new(1:L);
end
    
% match phase of output to input
if size(szy,2) == 3
    for latindx = 1:size(y,2)
        for lonindx = 1:size(y,1)
            
            IFinal = squeeze(Ymag(lonindx,latindx,:)).* sin(PhiX);
            RFinal = squeeze(Ymag(lonindx,latindx,:)).* cos(PhiX);
            
            spec = RFinal + 1i*IFinal;
            Y_NEW(lonindx,latindx,:) = spec;
            y_new(lonindx,latindx,:) = ifft(spec,'symmetric');
            
        end
    end
    
else
    if size(Ymag,2) == size(PhiX,1)
        Ymag = Ymag'; % fix so you dont get a matrix
        
    end
    IFinal = Ymag .* sin(PhiX);
    RFinal = Ymag .* cos(PhiX);
    spec= RFinal + 1i*IFinal;
    Y_NEW = spec;
    y_new = ifft(spec,'symmetric');
    
end

end

