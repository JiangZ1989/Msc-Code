function [S, F, T] = gen_spectrogram(x, window, noverlap, nsample, fs)
 N = length(x);
 S = [];
 pos = 1;
 while (pos+nsample <= N)
     frame = x(pos:pos+nsample-1);
     pos = pos + (nsample - noverlap);
     Y = fft(frame.*window, nsample); 
     S = [S Y(1:round(nsample/2), 1)]; %  mirror
 end
 [M, K] = size(S);
 F = (0:round(nsample/2)-1)' / nsample * fs; % [0, fs/2) Hz 
%   F = psdfreqvec(nsample, fs, 'half');
 T = (round(nsample/2):(nsample-noverlap):N-1-round(nsample/2))/fs;