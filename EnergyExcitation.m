%% Conner Brown
%  Date Created: 10/05/2016
%  Last Edited:  10/16/2016
%  File:    EnergyExcitation.m
%  Brief:   Perform frequency analysis on audio signal to determine bpm.
%           Use FFT on a small time interval and track when a large
%           spectral flux occurs. These locations in time correspond to
%           beats, from which, multiple may be extrapolated to bpm.

clear
% clc
format long g
directory=dir('Samples\*.mp3');
[s,~]=size(directory);
EEtest=cell(2,s);
% quarter notes at 159 bpm occur every .377 seconds or 16641 Samples
% FFT resolution must be precise enough to find these 
buf=2048;               % sample buffer for analyzing "real time"
fhistory=cell(5,1);     % subband average buffer
thistory=cell(5,1);     % time average buffer
bestbpm=0;              % what we are looking for

spec=zeros(buf,1);     % preparing spectrum array
mag=zeros(buf,1);      % preparing magnitude 
space=round(logspace(0,log(buf/2)/log(10),16));    % logarithmically spaced subbands
sub=zeros(15,1);        % preparing subband array
subsum=zeros(15,1);     % sum of subbands in history, for taking average and comparing to sub
factor=0;
beat=0;
prevbeat=0;
x=1;
reset=false;
slide=0;

% for 
file=32; % :s
    [y,Fs]=audioread(strcat('Samples\',directory(file).name));
    L=length(y);            % total number of samples in file
    t=(0:L-1)/Fs;           % time vecotr
    f=Fs*(0:(L/2))/L;       % frequency vector for whole length
    res=Fs/buf;             % resolution of frequency domain   
    for i=1:5
        fhistory{i,1}=zeros(15,1);
        thistory{i,1}=0;
    end
numbuf=(L-buf-mod(L,buf))/buf;
fullsum=0;                   % value of full spectrum energy
fullsumarr=zeros(numbuf,1);      % store full energy spectrum energy 
fullhistory=0;               % value of average full spectrum history energy
fullhistoryarr=zeros(numbuf,1);  % store average full spectrum history energy
timeenergy=0;                % value of time energy
timeenergyarr=zeros(numbuf,1);   % store time energy
tsum=1-1;               % value of average history time energy
thistoryarr=zeros(numbuf,1);  % store average time history energy, is average of thistory cell
peakdetect=zeros(3,1);          % compares first and third to second
lastpeakloc=0;          % Time energy location of old peak
lastpeakval=0;          % Time energy value of old peak
peakloc=0;              % Time energy location of new peak
peakval=0;              % Time energy value of new peak
tbpm=0;                 % Time Energy bpm every time we have a new beat
tbpms=zeros(1,1);
tbestbpm=0;             % this is what we are looking for Time Energy
tbpmindex=0;            % updates every window
subd=3*Fs;                 % smallest subdivision thus far
oldsubd=3*Fs;              % last subdivision to convert counts of integer multiples
intmult=0;              % integer multiple of smallest subdivision to get current subdivision
minreg=0;               % index of myints registered for current beat
minint=0;               % value of smallest difference between integer multiple and library of integer multiples
myints=(1:12)/12;   % array of empircally found integer multiples corresponding to subdivisions
mydiff=zeros(size(myints));              % array of difference between myints and intmult
myreg=zeros(size(myints));      % integer multiple register, keeps track of how often an integer multiple is found


tic
for window=1:buf:numbuf*buf
    reset=false;
    % Time Energy
    timeenergy=sum(y(window:window+buf-1,1).^2);
    timeenergyarr((window+buf-1)/buf,1)=timeenergy;
    % update history buffer with time energy value
    for i=1:4
        thistory{i,1}=thistory{i+1,1};
    end
    thistory{5,1}=timeenergy;
    for i=1:5
        tsum=((i-1)*tsum+thistory{i,1})/i;
    end
    thistoryarr((window+buf-1)/buf,1)=tsum;
    
    % peak detection
    for i=1:2
        peakdetect(i,1)=peakdetect(i+1,1);
    end
    peakdetect(3,1)=tsum;
    % middle must be greater than bookkends and greater than the
    % average...something how bout 85% of last peak value
    if peakdetect(2,1)>peakdetect(1,1) && peakdetect(2,1)>peakdetect(3,1) && peakdetect(2,1)>0.85*lastpeakval
        lastpeakloc=peakloc;
        lastpeakval=peakval;
        peakloc=window;
        peakval=tsum;
    end
    if peakloc==window
fprintf('Time excited!%d,%d\n',lastpeakloc,peakloc)
    end
    % keep track of smallest and next smallest subdivision
    if peakloc-lastpeakloc<subd
        oldsubd=subd;
        subd=peakloc-lastpeakloc;
        reset=true;
    end
    % if new subdivision then convert all registers
    if reset && oldsub<3*Fs
        slide=round(mod(oldsubd/subd*12,12));
        
        
    end
    if oldsub<3*Fs
        % see if current subdivision is an integer multiple of the subd
        intmult=(peakloc-lastpeakloc)/subd;
        % only look for meaningful subdivision integer multiples, determined by
        % time signature
        % keep track of how often the integer multiple is registered
        mydiff=myints-intmult;
        [minint,minreg]=min(abs(mydiff));
        % if minint is smaller than 2048 sampling error then register it
        if mydiff(minreg)>0
            if (peakloc-lastpeakloc+2048)/subd>myints(minreg)
                myreg(minreg)=myreg(minreg)+1;            
            end
        else
            if (peakloc-lastpeakloc-2048)/subd<myints(minreg)
                myreg(minreg)=myreg(minreg)+1;
            end
        end
    
        % convert sample distance to bpm
        [~,maxreg]=max(myreg);
        % tells me the most often played integer multiple of the
        % smallest subdivision, can extract time signature from this, as
        % well as true bpm
        myints(maxreg)*12;
        % average all bpms
        tbestbpm=((tbpmindex-1)*tbestbpm+tbpm)/tbpmindex;
        tbpmindex=tbpmindex+1;
    end
    % Frequency Energy
    spec=fft(y(window:window+buf-1,1));
    mag=abs(spec);
    %
    % why no phase? look into this...
    % 
    % concentrate magnitude spectrum into subbands using parsevals theorem
    for i=1:15
        sub(i,1)=(1/buf)*sum(mag(space(1,i):space(1,i+1),1).^2);
    end
    % full energy spectrum energy
    fullsum=sum(sub);
    fullsumarr((window+buf-1)/buf,1)=fullsum;
% %         %plot some magnitude spectrums
% %         figure
% %         plot(space(1,1:15),sub)
% %         title(window)
    % determine if this energy is much higher than average history energy
    subsum=zeros(15,1);
    for i=1:5
        subsum=((i-1)*subsum+fhistory{i,1})/i;
    end
    fullhistory=sum(subsum);
    fullhistoryarr((window+buf-1)/buf,1)=fullhistory;
    % threshold paramter, determined empirically
    if sub>1.3*subsum
        prevbeat=beat;
        beat=window;
        bpm=Fs*60/(beat-prevbeat);
        if bpm<80
           factor=ceil(80/bpm);
           if factor<10
                bpm=factor*bpm;
           end
        end
        if bpm>159
            factor=floor(bpm/80);
            if factor<10
                bpm=bpm/factor;
            end
        end
        % only look at the first excitement, delay until fastest quarter
        % note could happen
        % window=window+buf*8;
fprintf('Freq excited!%d,%d\n',prevbeat,beat)
    end
    % update history buffer with new subband array
    for i=1:4
    fhistory{i,1}=fhistory{i+1,1};
    end
    fhistory{5,1}=sub;
    % average all bpms
    bestbpm=((x-1)*bestbpm+bpm)/x;
    x=x+1;
end

toc

% plot time energy and frequency energies and averages
figure
hold on
plot(1:numbuf,timeenergyarr, 'r')
plot(1:numbuf,thistoryarr,'g')
plot(1:numbuf,fullsumarr, 'b')
plot(1:numbuf,fullhistoryarr,'y')
title(strcat(directory(file).name, ' excitations'))
legend('timeE', 'timeH','freqE','freqH')
hold off
fprintf('%s %d\n',directory(file).name,bpm)
% end % looping through whole directory