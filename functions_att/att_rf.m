
function [rf2] = att_rf(rf,x,z,fs)

[L1,L2] = size(rf)
c =1540;
for i = 1:L2
   line = rf(:,i); 
   
    t = z*2/c;
     figure; plot(t*1e6,line); xlabel('[\mus]')
     axis tight
     
     % x is signal; t is time vector.
% h = hilbert(line);
% unrolled_phase = unwrap(angle(h));
% inst_freq = diff(unrolled_phase)./diff(t)/(2*pi);
% inst_amp = abs(h);
%   
%figure; plot(line) 
%figure; plot(inst_freq)  
%figure; plot(inst_amp)

%y = exp(-c*t);
%figure


end



end