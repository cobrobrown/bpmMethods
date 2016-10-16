%% Conner Brown
%  Date Created 7/.../16
%  Last Edited  10/04/16
% BPM.m
% Find prominent "beats" in a song. Can be identified by rapid changes in
% spectral density (FFT) or repetitive spikes in amplitude (impulse train)
% This file serves as a repository for all ideas and attempts. See
% Improviser for implementation.
clear
clc
tic

directory=dir('Samples\*.mp3');
% for 
    file=35; % :46
%     if directory(file).bytes>10000000
%         file=file+1;
%     end
[y,Fs]=audioread(strcat('Samples\',directory(file).name));        % y is comprised of two channels, Fs is the sampling rate
L=length(y);
% t=(0:L-1)/Fs;
% f=Fs*(0:(L/2))/L;

%% Plotting the waveform in time space
% tic
% figure
% % subplot(2,1,1)
% plot(t,y(:,1))
% title(strcat('Channel 1 of ',directory(1).name))
% xlabel('time (s)')

% subplot(2,1,2)
% plot(t,y(:,2))
% title('Channel 2')
% xlabel('time (s)')

% toc
%% Determining the bpm
% % Method 1:Test peak to peak of set relative to last set for rapid change in
% % amplitude. 
% tic
% lastind=0;
% count=1;
% index=0;
% space=0;
% bpm=0;
% testavbpm=0;
% avbpm=0;
% beginind=0;
% step=ceil(.03*Fs);
% for index=1:step:(L-2*step)      % Proceed in 0.01 seconds' worth of samples
%     [maxpk,maxind]=max(y(index:(index+step),1));
%     maxind=maxind+index;
%     [minpk,minind]=min(y(index:(index+step),1));
%     minind=minind+index;
%     [nextmaxpk,nextmaxind]=max(y((index+step):(index+2*step),1));
%     nextmaxind=nextmaxind+index+step;
%     [nextminpk,nextminind]=min(y((index+step):(index+2*step),1));
%     nextminind=nextminind+index+step;
%     pkpk=maxpk-minpk;
%     nextpkpk=nextmaxpk-nextminpk;
%     if nextpkpk>3*pkpk && nextpkpk>1.2      % criterion parameters...slightly arbitrary
%         if count>2
%             space=nextmaxind-beginind;
%             for multiple=1:4                % see if spacing is a multiple of an ordinary bpm
%                 bpm=space/Fs*60*multiple;
%                 if bpm>61 && bpm<179
%                     avbpm=(avbpm+bpm)/2;
%                     break
%                 end
%                 bpm=space/Fs*60/multiple;
%                 if bpm>61 && bpm<179
%                     avbpm=(avbpm+bpm)/2;
%                     break
%                 end
%             end
%         end
%         beginind=nextmaxind;
%         count=count+1;
%         index=index+Fs/4*180/60;
%     end
% %     if count>20
% %         break
% %     end
% end
% 
% toc
% % Method 2: look for peaks above a threshold, find position of peaks above threshold and find distances
% % between peaks which are contained in limits of normal bpms (61-179)
% tic
% lastind=0;
% count=1;
% index=0;
% space=0;
% bpm=0;
% avbpm=0;
% % Find the threshold value. Max bpm should be 180 so in a song x minutes
% % long there should only be 180*x beats total. Look for a reasonable number
% % for these peaks, the lowest value fo these is the threshold.
% minbeats=180/60*L/Fs;
% % Find the spacing (then bpm) of these threshold peaks
% for index=1:L
%     if y(index,1)>0.4
%         space=index-lastind;
%         lastind=index;
%         for multiple=1:4
%             bpm=space/Fs*60*multiple;
%             if bpm<179 && bpm>61
%                 avbpm=(avbpm+bpm)/2;
%                 count=count+1;
%                 break
%             end
%         end
%         index=index+Fs/4; %minimum bpm is 60 so skip everything in between
%     end
%     if count>40
%         break
%     end
% end
% 
% toc
% Method 3: Impulse Response correlation. Generate impulse responses at
%   bpms ranging from 80 to 159. Convolve each one with the time waveform,
%   our bpm is the one with highest correlation to impulse.
tic
high=0;                 % placeholder to compare every iteration of phase
bestbpm=0;              % holds value of what we're trying to find
best=0;                 % placeholder to compare every iteration of bpm
test=0;
sum=0;                  % sums value of waveform at impulse loaction
bpmwidth=0;             % sample distance between impulses
avgbpm=0;               % running average of bpms for window

for bpm=80:159
    % tic
    bpmwidth=ceil(Fs*60/bpm);
     % generate impulse train
%     train=zeros((3*bpmwidth),1);         
%     for value=1:bpmwidth:length(train)
%         train(value,1)=1;  
%     end
%     for i=1:bpmwidth
%         test=dot(train,y(((i+30*Fs):(i+length(train)-1+30*Fs)),1));       % Hopefully no boring ass intro
%         if test>high
%             high=test;
%         end
%     end
    %  instead, just pick out values where impulse train is
    for phase=1:ceil(bpmwidth/100):bpmwidth        
        % Change step size to larger values to make code run faster, seems 
        % bpmwidth/25 is upper limit, bpmwidth/100 is appropriately
        % accurate
        sum=0;
        % sum at impulse location
        for i=phase:bpmwidth:(phase+(L-2*bpmwidth))
            sum=sum+y(i,1);
        end
        % update highest sum for phase
        if abs(sum)>high
            high=abs(sum);
        end
    end
    % update highest sum for bpm
    if high>best
        best=high;
        bestbpm=bpm;
    end
end
% end
toc
% clearvars i j bpmwidth test phase % sum high bpm
fprintf('%s %d\n',directory(file).name,bestbpm)
% end    % For running through the whole directory 

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

fprintf('%s %d\n',directory(file).name,bestbpm)
toc


%% Short Time Fourier Transform
% Method 1: spectrogram, gives short time FT info specifically when the FT
% rapidly changes in short time spans
% tic
% spectrogram(y(:,1),'yaxis')
% 
% toc
% % 
% % % Method 2: localized FFT, same idea as Method 1 but I choose convenient
% % % locations to do the FFT
% % starttime=10;        % in seconds, should happen exactly on beat
% % around=.37;            % in seconds, how far away to look between 0.37 for 159 bpm and 0.75 for 80 bpm       
% % start=starttime*Fs; % convert to samples
% % sta=ceil(start-around*Fs/2);
% % fin=ceil(start+around*Fs/2);
% % % length must be between that allowed for bpms
% % % length<33283 for 159 bpm & length>66150 for 80 bpm
% % length=fin-sta;
% % freq=Fs*(0:(length/2))/length;
% % shFFT=fft(y(sta:fin,:));
% % fullSpec=abs(shFFT/length);
% % chanleft=fullSpec(1:ceil(length/2),1);
% % chanleft(2:end-1)=2*chanleft(2:end-1);
% % % plot that ish
% % figure
% % plot(freq,chanleft)
% % title('Short FFT')
% % 
% % toc


%% Differentiating waveform
% % % % load('tempdata.mat')
% % dydt=diff(y(:,1))/(L/Fs);     % Basically deltay/deltax
% % dydt(length(dydt)+1,1)=0;     % add a zero at the end
% % % ddydt=diff(dydt)/(L/Fs);
% % % ddydt(length(ddydt)+1,1)=0;
% % dy2=dydt.^2;
% % % y2=y(:,1).^2;
% % % ddy=diff(y2)/(L/Fs);
% % % ddy(length(ddy)+1,1)=0;
% % % figure 
% % % plot(t,y(:,1))
% % % figure 
% % % plot(t,dydt)
% % 
% % % Method 3: Impulse Response correlation. Generate impulse responses at
% % %   bpms ranging from 80 to 159. Convolve each one with the time waveform,
% % %   our bpm is the one with highest correlation to impulse.
% % tic
% % high=0;                 % placeholder to compare every iteration of phase
% % bestbpm=0;              % holds value of what we're trying to find
% % best=0;                 % placeholder to compare every iteration of bpm
% % test=0;
% % sum=0;                  % sums value of waveform at impulse loaction
% % bpmwidth=0;             % sample distance between impulses
% % for bpm=80:159
% %     % tic
% %     bpmwidth=ceil(Fs*60/bpm);
% %      
% % %     train=zeros((3*bpmwidth),1);         % generate impulse train
% % %     for value=1:bpmwidth:length(train)
% % %         train(value,1)=1;  
% % %     end
% % %     for i=1:bpmwidth
% % %         test=dot(train,y(((i+30*Fs):(i+length(train)-1+30*Fs)),1));       % Hopefully no boring ass intro
% % %         if test>high
% % %             high=test;
% % %         end
% % %     end
% %     for phase=1:ceil(bpmwidth/100):bpmwidth        % Change step size to larger vlues to make code run faster, seems bpmwidth/25 is limit
% %         sum=0;
% %         for i=phase:bpmwidth:(phase+(L-2*bpmwidth))
% %             sum=sum+dy2(i,1);
% %         end
% %         if abs(sum)>high
% %             high=abs(sum);
% %         end
% %     end
% %     if high>best
% %         best=high;
% %         bestbpm=bpm;
% %     end
% %     % toc
% % end
% % 
% % toc

% % %% FFT of waveforms
% % tic
% % % fft of 2D vector y
% % yfft=fft(y);
% % full=abs(yfft/L);
% % % channel 1
% % chan1=full(1:L/2+1,1);
% % chan1(2:end-1)=2*chan1(2:end-1);
% % % channel 2
% % % chan2=full(1:L/2+1,2);
% % % chan2(2:end-1)=2*chan2(2:end-1);
% % % plotting fft
% % figure 
% % % subplot(2,1,1)
% % plot(f,chan1)
% % title('channel1 spectrum')
% % % subplot(2,1,2)
% % % plot(f,chan2)
% % toc

%% Bass Inquiry

% % %% Method 1:
% % %  Filter out high pitches, look just for bass, recombine and determine bpm
% % tic
% % df=designfilt('lowpassfir','PassbandFrequency',100/Fs,'StopbandFrequency',800/Fs);
% % D=mean(grpdelay(df));
% % filtered=filter(df,y(:,1));
% % figure
% % plot(t,filtered)
% % title('filtered wave')
% % filtfft=fft(filtered);
% % filtfull=abs(filtfft/L);
% % filtspec=filtfull(1:L/2+1,1);
% % filtspec(2:end-1)=2*filtspec(2:end-1);
% % figure
% % plot(f,filtspec)
% % title('filter spectrum')
% % toc
% % 
% % % % Method 3: Impulse Response correlation. Generate impulse responses at
% % % %   bpms ranging from 80 to 159. Convolve each one with the time waveform,
% % % %   our bpm is the one with highest correlation to impulse.
% % % tic
% % % high=0;                 % placeholder to compare every iteration of phase
% % % bestbpm=0;              % holds value of what we're trying to find
% % % best=0;                 % placeholder to compare every iteration of bpm
% % % test=0;
% % % sum=0;                  % sums value of waveform at impulse loaction
% % % bpmwidth=0;             % sample distance between impulses
% % % for bpm=80:159
% % %     bpmwidth=ceil(Fs*60/bpm);
% % %     for phase=1:ceil(bpmwidth/100):bpmwidth        % Change step size to larger vlues to make code run faster, seems bpmwidth/25 is limit
% % %         sum=0;
% % %         for i=phase:bpmwidth:(phase+(L-2*bpmwidth))
% % %             sum=sum+filtered(i,1);
% % %         end
% % %         if abs(sum)>high
% % %             high=abs(sum);
% % %         end
% % %     end
% % %     if high>best
% % %         best=high;
% % %         bestbpm=bpm;
% % %     end
% % % end
% % % 
% % % toc

% %% Method 2:
% % %  Inverse FFT truncating high frequencies (>100 Hz) i.e. lower the sample
% % %  rate of the time waveform to 100 Hz
% % 
% % % FFT of waveforms
% % tic
% % % fft of 2D vector y
% % yfft=fft(y,L);
% % full=abs(yfft/L);
% % % channel 1
% % chan1=full(1:L/2+1,1);
% % chan1(2:end-1)=2*chan1(2:end-1);
% % % channel 2
% % % chan2=full(1:L/2+1,2);
% % % chan2(2:end-1)=2*chan2(2:end-1);
% % % plotting fft
% % figure 
% % % % subplot(2,1,1)
% % plot(f,chan1)
% % title('channel1 spectrum')
% % % subplot(2,1,2)
% % % plot(f,chan2)
% % toc
% % 
% % 
% % % Not necessary, literally jut make everything above 100 Hz 0...
% % % Take <100 Hz
% % % y100=yfft((ceil(length(y(:,1))/2)-100):(ceil(length(y(:,1))/2)+100),1);
% % % invy=ifft(y100,L,'symmetric');
% % % figure
% % % plot(t,invy)
% % % title('undersampled')
% % rudefilt=chan1;
% % rudefilt(ceil((100/Fs*L)):end,1)=0;
% % % figure
% % % plot(f,rudefilt)
% % invfilt=ifft(rudefilt,L,'symmetric');
% % figure
% % plot(t,invfilt)
% % 
% % % Method 3: Impulse Response correlation. Generate impulse responses at
% % %   bpms ranging from 80 to 159. Convolve each one with the time waveform,
% % %   our bpm is the one with highest correlation to impulse.
% % tic
% % high=0;                 % placeholder to compare every iteration of phase
% % bestbpm=0;              % holds value of what we're trying to find
% % best=0;                 % placeholder to compare every iteration of bpm
% % test=0;
% % sum=0;                  % sums value of waveform at impulse loaction
% % bpmwidth=0;             % sample distance between impulses
% % for bpm=80:159
% %     bpmwidth=ceil(Fs*60/bpm);
% %     for phase=1:ceil(bpmwidth/100):bpmwidth        % Change step size to larger vlues to make code run faster, seems bpmwidth/25 is limit
% %         sum=0;
% %         for i=4*phase:bpmwidth:(phase+(L-4*bpmwidth))
% %             sum=sum+invfilt(i,1);
% %         end
% %         if abs(sum)>high
% %             high=abs(sum);
% %         end
% %     end
% %     if high>best
% %         best=high;
% %         bestbpm=bpm;
% %     end
% % end
% % 
% % toc
% % % looking at periodicity in filtered time wave 
% % avgfilt=mean(invfilt(Fs:length(invfilt)-Fs,1).^2);

%% Prize-winning BPM
%  Look deeply into FFT of waveform (0-30 Hz). Peaks indicate rhythms and
%  tempos. BPMs between 80-159 correspond to 1.33-2.65 Hz. Triple, quarter
%  and eighth notes will be at 3,4,and 8 multiples of BPM. So, Looking at
%  3*1.33-8*2.65 (or 4.0-21.2 Hz). These peaks give BPM and Meter.

% FFT of waveforms
tic
% fft of 2D vector y
yfft=fft(y);
full=abs(yfft/L);
clear yfft
% channel 1
chan1=full(1:L/2+1,1);
clear full
chan1(2:end-1)=20*log10(2*chan1(2:end-1));
% plotting fft
% figure 
% plot(f,chan1)
% title('channel1 spectrum')
% axis([0 30 0 3*10^-3])
% Magnitude in dB
% plot(f,chan1)
% title('Chan1 Mag (dB)')

% find peaks
% Method 1: Look at fine resolution sweep
% % [pks,locs]=findpeaks(chan1(1:ceil(30*L/Fs),1));
% % % only look for localized peaks (1 BPM = 0.0167 Hz)
% % localpk=zeros(1270,2);
% % for index=floor(((4/3-0.167/2):0.0167:(8*159/60-0.0167/2))*L/Fs);
% %     [localpk(index,1),localpk(index,2)]=max(chan1(index:index+(0.0167/2)*L/Fs,1));
% % end
% % % compare local peaks to find prominent meter/bpm

% Method 2: Look at bpm intervals
% look at n*1.33:2.65 intervals. peaks divided by interval number (meter)
% is equal to bpm
firstbpm=round(80/60*L/Fs);
intervalpks=zeros(8,2);
bpms=zeros(8,1);
for interval=1:8
    [intervalpks(interval,1),intervalpks(interval,2)]=max(chan1(interval*firstbpm:(interval+1)*firstbpm,1));
    bpms(interval,1)=(intervalpks(interval,2)+interval*firstbpm-1)/interval*Fs/L*60;
end
trustBPM=bpms(1,1);
wholeBPM=mean(bpms);
wholestd=std(bpms);
% weighting of FFT values can be used from intervalpk(:,1)

% Do this every 5 seconds to see if certain parts of the song suppress BPM
tic
sec=Fs*60;
smallfirstbpm=round(4/3*sec/Fs);
fftpks=zeros(8,2);
smallbpms=zeros(8,ceil((L-sec)/sec));
for n=1:sec:(L-sec)
    smallfft=fft(y(n:n+sec,1));
    smallfull=abs(smallfft/sec);
    clear smallfft
    smallspec=smallfull(1:sec/2+1,1);
    clear smallfull
    smallspec(2:end-1)=2*smallspec(2:end-1);
    for interval=1:8
        [fftpks(interval,1),fftpks(interval,2)]=max(smallspec(interval*smallfirstbpm:(interval+1)*smallfirstbpm,1));
        smallbpms(interval,ceil(n/sec))=(fftpks(interval,2)+interval*smallfirstbpm-1)/interval*Fs/sec*60;
    end
end
toc
% figure
% smallfreq=(0:(sec/2))*Fs/sec;
% plot(smallfreq,smallspec)
% title('small spectrum')
A=mean(smallbpms);
smallBPM=mean(A);
smallstd=std(A);

fprintf('%s %d %d %d %d\n',directory(file).name,bestbpm,wholeBPM,trustBPM,smallBPM)
% end