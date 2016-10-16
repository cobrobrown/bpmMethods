%% Conner Brown
%  Date Created: 10/04/2016
%  Last Edited:  10/05/2016
%  File:    CombFilter.m
%  Brief:   Use impulse train to simulate a digital comb filter in audio 
%           data. Cycle through select impulse trains representing beats at
%           certian bpm. After passing audio data through comb filter, 
%           audio bpm should emerge as highest sum of waveform values.

% Impulse Train, Digital Comb Filter.
%           The idea is to generate an impulse train with spacing
%           determined by bpm (60 bpm = 1 impulse per 44100 Samples). We 
%           take the dot product of the impulse train and the audio data to
%           test how strongly correlated they are. We shift the impulse 
%           train one sample and take the dot product again. Do this until 
%           we've covered the spacing of this impulse train. Whichever dot 
%           product is highest, we store as this bpm's highest value. Do 
%           this for every bpm. The highest dot product value of all bpms 
%           is the winner! 
%           This models a digital comb filter by only looking at values
%           located at the impulses and attenuating (setting to 0) 
%           everything else.
%           My code below does not generate an impulse train and take the
%           dot product as this would be computationally hefty. Instead, I
%           sum the audio data values which the impulse train would have
%           picked out from doing the dot product. Along the lines of
%           saving on computation time, I also shift the "impulse train"
%           more than 1 sample at a time. Empirical values are in comments.
%           Since my testing set consists of full songs, but the
%           application in mind works with realtime data input, I work with
%           5 second increments. A running average is computed from the
%           start of the song as if that is when the program went live.

clear
% clc
format long g
directory=dir('Samples\*.mp3');
[s,~]=size(directory);
testarray=cell(2,s);
% for 
    file=33; % :s
[y,Fs]=audioread(strcat('Samples\',directory(file).name));        % y is comprised of two channels, Fs is the sampling rate
L=length(y);                % Length of audio file in Samples
t=(0:L-1)/Fs;               % time vector
f=Fs*(0:(L/2))/L;           % frequency vector

tic
high=0;                 % placeholder to compare every iteration of phase
bestbpm=0;              % holds value of what we're trying to find
best=0;                 % placeholder to compare every iteration of bpm
test=0;
sum=0;                  % sums value of waveform at impulse loaction
bpmwidth=0;             % sample distance between impulses
tic
% starting point for analysis on y, 3 second update
for window=1:(3*44100):(L-3*44100-mod(L,3*44100));
    for bpm=80:159
        bpmwidth=ceil(Fs*60/bpm);
        quarternote=ceil(Fs*60/bpm/4);
        for phase=1:ceil(quarternote/25):quarternote  
            quartersum=0;
            % in 3 seconds, we can see at most 15 quarter notes at 80 bpm
            for i=1:15
                quartersum=quartersum+y(window+i*quarternote+phase,1);
            end
            if abs(quartersum)>high
                high=abs(quartersum);
            end
        end 
        if high>best
            best=high;
            bestbpm=bpm;
        end       
    end
end
time=toc;

clearvars i j bpmwidth test phase % sum high bpm
fprintf('%s %d\n %d\n',directory(file).name,bestbpm,time)
% temp=directory(file).name;
% temp=str2num(temp(1:2));
% if temp<60
%   temp=temp+100;
% end
% testarray{1,file}=temp;
% testarray{2,file}=bestbpm;
% end    % For running through the whole directory 

%% Proving Parseval's Theorem
my=y(1:44100,1);
timesum=0;
for i=1:length(my)
    timesum=timesum+my(i)^2;
end
freqsum=0;
% 1024 samples = 0.001436 seconds
% 44100 samples = 0.003667 seconds
tic
spec=fft(my);
mymag=abs(spec);
toc
% myfft=mymag(1:length(my)/2+1,1);
% myfft(2:end-1)=2*myfft(2:end-1);
for i=1:length(myfft)
    freqsum=freqsum+myfft(i)^2;
end