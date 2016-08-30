function [S, F, T] = spStft(x, window, frame_overlap, frame_length, fs)
 nsample = round(frame_length  * fs / 512); % convert ms to points
 noverlap = round(frame_overlap * fs / 512); % convert ms to points
 window   = eval(sprintf('%s(nsample)', window)); % e.g., hamming(nsample)
 %% spectrogram
 [S, F, T] = gen_spectrogram(x, window, noverlap, nsample, fs); % below
 %% plot
%  if show
%      plot_spectrogram(S, F, T); % below
%  end