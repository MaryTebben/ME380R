%% This set of files was found on 
% https://www.clear.rice.edu/elec301/Projects01/beat_sync/beatalgo.html

% CONTROL takes in the names of two .wav files, and outputs their
% combination, beat-matched, and phase aligned.
%
%     SIGNAL = CONTROL(SONG1, SONG2, BANDLIMITS, MAXFREQ) takes in
%     the names of two .wav files, as strings, and outputs their
%     sum. BANDLIMITS and MAXFREQ are used to divide the signal for
%     beat-matching
%
%     Defaults are:
%        BANDLIMITS = [0 200 400 800 1600 3200]
%        MAXFREQ = 4096
  
%   if nargin < 3, bandlimits = [0 200 400 800 1600 3200]; end
%   if nargin < 4, maxfreq = 4096; end
%   song2 = ('04.DavidBowieQueen-UnderPressure.mp3');
  clc; clear; close all;
  format compact;

  % Input desired song as an mp3 file
  song1 = ('08.Queen-WeWillRockYou-weAreTheChampions.mp3');  

  % audioread finds the audio point information (x1) and the song frequency
  % (sr1)
  [x1,sr1] = audioread(song1);
  song = x1;

  % maxfreq gives the maximum band limit for audio type separation
  maxfreq = sr1/2;
  sample_size = floor(2.2*2*maxfreq); 
  bandlimits = [0 200 400 800 1600 3200];
  
  start = sr1*23.75; % start time for audio clip
  stop = sr1*29.6;   % end time for audio clip
  
  sample = song(start:stop);
  figure(5)
  plot(sample)
  dss = fft(sample);
  dss_re = abs(dss);
  figure(6)
  plot(dss_re)
  % Implements beat detection algorithm for each song
  
  %% implementation given by Rice University
  status = 'filtering first song...'
  a = filterbank(sample, bandlimits, maxfreq);
  status = 'windowing first song...'
  b = hwindow(a, 0.2, bandlimits, maxfreq);
  status = 'differentiating first song...'
  c = diffrect(b, length(bandlimits));
  status = 'comb filtering first song...'
  %    Recursively calls timecomb to decrease computational time
  d = timecomb(c, 2, 60, 240, bandlimits, maxfreq);
  e = timecomb(c, .5, d-2, d+2, bandlimits, maxfreq);
  f = timecomb(c, .1, e-.5, e+.5, bandlimits, maxfreq);
  g = timecomb(c, .01, f-.1, f+.1, bandlimits, maxfreq);
  song_bpm = g;

  figure(1)
  plot([1:length(sample)], abs(a))
  %%

  b_simp = b(:,1)+b(:,2)+b(:,3)+b(:,4)+b(:,5);

  figure(2)
  plot([1:length(sample)], abs(b))
  legend();
  figure(3)
  plot([1:length(sample)], abs(b(:,2)))
  hold on
  % finds peaks in desired audio channel
  findpeaks(b(:,2),1:length(b(:,2)),'MinPeakProminence',4,'Annotate','extents');
  findpeaks(b(:,2),1:length(b(:,2)),'MinPeakDistance',10000);
  % minPeakDistance of 10,000 because 44100 points is 1 second
%   hold off

  % multiple peak locatio methods for comparison purposes
  [b_peaks,b_loc] = findpeaks(b(:,2),1:length(b(:,2)),'MinPeakProminence',4,'Annotate','extents');
  [b_peaks_sr,b_loc_sr] = findpeaks(b(:,2),1:length(b(:,2)),'MinPeakDistance',10000);
  peak = [b_peaks b_loc'];
  peak_sr = [b_peaks_sr b_loc_sr'];

  peaks = intersect(peak,peak_sr,'rows');
  peak = sortrows(peaks,2);

  j = 1;
  k = 1;
  l = 1;
  m = 1;
  breakpoint = length(b(:,2))/4;
  % searches through each defined section of the song and removes obviously
  % incorrect peaks
  for i = 1:length(peak)
      if peak(i,2) < breakpoint && peak(i,1) > 200
          peaks1(j,:) = peak(i,:);
          j = j+1;
      elseif peak(i,2) > breakpoint && peak(i,2) < 2*breakpoint && peak(i,1) > 175
          peaks2(k,:) = peak(i,:);
          k = k+1;
      elseif peak(i,2) > 2*breakpoint && peak(i,2) < 3*breakpoint && peak(i,1) > 150
          peaks3(l,:) = peak(i,:);
          l = l+1;
      elseif peak(i,2) > 3*breakpoint && peak(i,2) < 4*breakpoint && peak(i,1) > 90
          peaks4(m,:) = peak(i,:);
          m = m+1;
      end
  end

  peak = [peaks1; peaks2; peaks3; peaks4];
  peaks = intersect(peak,peak_sr,'rows');
  peaks = sortrows(peaks,2);

  plot(peaks(:,2),peaks(:,1), '*', 'linewidth', 3)
%   hold off

%% Take found peaks and turn them into desired audio
% entirely tunable

  % removes peaks that would interfere with desired tune (even if they are
  % correct peak locations)
  edited(:,:) = peaks(:,:);
  for i = 1:length(edited)
      if i == 2 || i == 4
          edited(i,:) = [];
      end
  end
    plot(edited(:,2),edited(:,1), 'o', 'linewidth', 3)
    hold off

  edited = sortrows(edited,2);
  cam = zeros(length(b),1);
  cam = [cam (1:length(b(:,1)))'];

  j = 1;
  k = 1;
  l = 1;
  m = 1;
  % separates peak locations into defined sections of the song
  for i = 1:length(edited)
    if edited(i,2) < breakpoint 
        edited1(j,:) = edited(i,:);
        edited1(j,1) = -1;
        j = j+1;
    elseif edited(i,2) > breakpoint && edited(i,2) < 2*breakpoint
        edited2(k,:) = edited(i,:);
        edited2(k,1) = 1;
        k = k+1;
    elseif edited(i,2) > 2*breakpoint && edited(i,2) < 3*breakpoint
        edited3(l,:) = edited(i,:);
        l = l+1;
    else %edited(i,2) > 3*breakpoint && edited(i,2) < 4*breakpoint
        edited4(m,:) = edited(i,:);
        m = m+1;
    end
  end

  % the minimum peak in each section should be a different note to the
  % reminaing peaks
  minimum3 = min(edited3);
  minimum4 = min(edited4);

  % -1 and 1 determine whether the note correlates to a protrusion or
  % cavity in the cam
  for i = 1:length(edited3)
    if edited3(i,1) > minimum3(:,1)
        edited3(i,1) = 1;
    else
        edited3(i,1) = -1;
    end
  end
  for i = 1:length(edited4)
    if edited4(i,1) > minimum4(:,1)
        edited4(i,1) = 1;
    else
        edited4(i,1) = -1;
    end
  end

  edited = [edited1;edited2;edited3;edited4];
  for i = 1:length(edited)
      j = edited(i,2);
      cam(j,:) = edited(i,:);
  end
%   cam = sortrows(cam,2);
  cam = [cam(:,1) cam(:,2)/44100];
  figure(4)
  plot(cam(:,2),cam(:,1),'-')

  circle_rad = 70; %mm
  circum = pi*2*circle_rad;

  % locates audio peaks around desired cam nominal circumference 
  max_time = cam(length(cam),2);
  scale = circum/max_time;
  cam = [cam(:,1) cam(:,2)*scale];
  angles = linspace(0,2*pi,circum);

  x = circle_rad*cos(angles);
  y = circle_rad*sin(angles);

  figure(7)
  axis square; grid on;
  plot(x,y)
  hold on 
  nonzero = find(cam(:,1));
  nonzero_angle = nonzero*scale/44100*(2*pi/circum);
  x_nonzero = circle_rad*cos(nonzero_angle);
  y_nonzero = circle_rad*sin(nonzero_angle);
  for i = 1:length(edited)
      if edited(i,1) == -1
        plot(x_nonzero(i),y_nonzero(i),'ko','linewidth',2)
      elseif edited(i,1) == 1
        plot(x_nonzero(i),y_nonzero(i),'ro','linewidth',2)
      end
  end
  hold off

  % finds angle between peaks for solidworks cam-creation purposes
  complex = x_nonzero+1i*y_nonzero;
  angle_comp = rad2deg(angle(complex))

