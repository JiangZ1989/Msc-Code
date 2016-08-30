clear all;
close all;
% 
% %Connecting to Outlook
% outlook = actxserver('Outlook.Application');
% mapi=outlook.GetNamespace('mapi');
% INBOX=mapi.GetDefaultFolder(6);
% 
% % Retrieving last email
% count = INBOX.Items.Count; %index of the most recent email.
% firstemail=INBOX.Items.Item(count); %imports the most recent email
% % secondmail=INBOX.Items.Item(count-1); %imports the 2nd most recent email
% subject = firstemail.get('Subject');
% body = firstemail.get('Body');
% 
% % Saving attachments to current directory
% attachments = firstemail.get('Attachments');
% if attachments.Count >=1
%     fname = attachments.Item(1).Filename;
%     dir = pwd;
%     full = [pwd,'\',fname];
%     attachments.Item(1).SaveAsFile(full)
% end


[x, fs] = audioread('object1.wav');
t=(0:length(x)-1)/fs; %time axis
L=length(x);
x = x + 0.01*rand(size(x)); 
% sigma=0.02;  %%could be change
% noise = x + sigma*randn(size(x));
% x=noise+x;
figure(1);
subplot(4,1,1);
plot(t,x);
xlabel('time(sec)'); %plot time domian
figure(1);
subplot(4,1,2);

NFFT = 2^(nextpow2(L));
Yfft=fft(x,NFFT)/L;
f = fs/2*linspace(0,1,NFFT/2+1);
plot(f,2*abs(Yfft(1:NFFT/2+1))); 
delta_f=fs/NFFT;

%  ydb=20*log10(y);
% xlabel('Frequency(Hz)');ylabel('magnitude(dB)');title('Frequency domain');
xlabel('Frequency(Hz)');ylabel('amplitude');title('Frequency domain');

%apply kaiser bandpass filter
[n,wn,beta] = kaiserord([50 100 4000 4200],[0 1 0],[0.05 0.001 0.05],44100);
h_wind = fir1(n,wn,kaiser(n+1,beta),'noscale');% impluse response coeffienct 
%fvtool(h_wind,1);%freqz(h_wind,1); 

global buffer;
buffer = zeros(size(h_wind));
for ii=(1:length(x))
 output(ii)=FIR(x(ii),h_wind);
end

subplot(4,1,3);
plot(output);
title('filted sounds');
xlabel('time/s');
Y=output;
 
Y_mov = tsmovavg(Y, 's', 1000, 2); 
subplot(4,1,4);
plot(Y_mov);
title('moving average');
xlabel('time(s)');
 
% fs = 44100;
frame_overlap = 10; % ms
frame_length  = 20;
figure(2);
%STFT
[S, F, T]=spStft(Y_mov(:), 'hanning', frame_overlap, frame_length, fs);
SdB = 20 * log10(abs(S)); % dB

%============The above code is from Yingying Tang==============

A=S.';
P=real(A).^2 + imag(A).^2; %power spectrum
[LT LF]= size(P); %Time and frequency segement numbers;

% find peaks on power spectrum
p = 1;
for t=1:LT
    for f = 2:LF-1
        if (P(t,f)>P(t,f-1) && P(t,f)>P(t,f+1))
            Ppeak(t,p) = P(t,f);
            kpeak(t,p) = f;
            p = p+1;
        end
    end
    Np(t) = p-1;
    p=1;
end

Bcrest = 8;%band of spectrum subsets

%calculate the mean and standard deviation of each subset
for t = 1:LT
    for p = 1:Np(t)
        if ((kpeak(t,p)>(Bcrest/2)) && (kpeak(t,p)+Bcrest/2<=LF))
            MeanPband(t,p) = mean(P(t,kpeak(t,p)-Bcrest/2 : kpeak(t,p)+Bcrest/2));
            StdevPband(t,p) = std(P(t,kpeak(t,p)-Bcrest/2 : kpeak(t,p)+Bcrest/2));
        else
            MeanPband(t,p) = 65536;
            StdevPband(t,p) = 65536;
        end
    end
end
 
%constants
Cmean = 0.9;
Cstd = 0.9;

% find crests on the power spectrum
for t = 1:LT
    c = 1;
    for p = 1:Np(t)
        if Ppeak(t,p) > Cmean*MeanPband(t,p) + Cstd*StdevPband(t,p)
            Pcrest(t,c) = Ppeak(t,p);
            kcrest(t,c) = kpeak(t,p);
            c= c+1;
        end
    end
    Nc(t) = c-1; 
end

%crests' indicator
B = zeros(size(P));
for t = 1:LT
    for c = 1:Nc(t)
        if kcrest(t,c) > 0
            B(t,kcrest(t,c)) = 1; %1 indicates peaks and 0 indicates no 
        end
    end
end

figure(2);
imshow(flipud(B'));
xlabel('Time');
ylabel('Frequency');
title('Crests that potentially originated by wheezing');
    

C=B;
delta_t=0.01;
kmin=0.1/delta_t;  % 100ms
k=2;
for i=2:LT-1  %%lable the same wheeze
    for j=2:LF-1
        if C(i-1,j+1)==0 && C(i-1,j)==0 && C(i-1,j-1)==0 && C(i,j)==1
            C(i,j)=k;  %% lable the new start of a wheeze
            k=k+1;
        end
        if (C(i,j)>1)
            if C(i+1,j+1)==1 
                C(i+1,j+1)=C(i,j);
            end
                if C(i+1,j)==1
                    C(i+1,j)=C(i,j);
                end
                    if C(i+1,j-1)==1
                        C(i+1,j-1)=C(i,j);
                    end
        end
     end
end

for n=1:max(C(:))    %% evaluate the length of the wheeze (100ms) 
    [row,col]=find(C==n);
    length_time=max(col)-min(col)+1;
    if length_time < 10
        C(C==n)=0;
    end
end

if max(C(:))==0    %% display whether it detected wheeze
    disp('no wheeze detected')
else
    disp('detected wheeze in the sample')    
end

%Display wheezes
figure(3);
imshow(flipud(C'));
xlabel('Time');
ylabel('Frequency');
title('Detected Wheezes');



% 
% %send test result to user
% mail = 'tangyingyingyy@gmail.com'; %Your GMail email address
% password = 'ping@070709';  %Your GMail password
% setpref('Internet','SMTP_Server','smtp.gmail.com');
% setpref('Internet','E_mail',mail);
% setpref('Internet','SMTP_Username',mail);
% setpref('Internet','SMTP_Password',password);
% props = java.lang.System.getProperties;
% props.setProperty('mail.smtp.auth','true');
% props.setProperty('mail.smtp.socketFactory.class', 'javax.net.ssl.SSLSocketFactory');
% props.setProperty('mail.smtp.socketFactory.port','465');
% 
% % add= 'guigu1289@gmail.com' 
% if max(C(:))==0    %% display whether it detected wheeze  
%     sendmail(useraddress, 'test result', 'no wheeze detected.','results.fig');
% else  
%     sendmail(useraddress, 'test result', 'detected wheeze in the sample','results.fig');
% end
% 

