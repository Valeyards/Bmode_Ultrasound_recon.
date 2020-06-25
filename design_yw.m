%% clear workspace
clc;clear all;
%% load_data
 img1 = LoadRF('data/RF_B_00000.pgm');
img2 = LoadRF('data/RF_B_00001.pgm');
img3 = LoadRF('data/RF_B_00002.pgm');
img4 = LoadRF('data/RF_B_00003.pgm');
img5 = LoadRF('data/RF_B_00004.pgm');
img6 = LoadRF('data/RF_B_00005.pgm');
f=(img1+img2+img3+img4+img5+img6)/6;
% f1=f';
img=f;
% figure,imagesc(abs(f)),colormap('gray');colorbar;
% %[num_lines,mgraph]=size(f);

%% 参数设置
[lines,points]=size(img);
depth=3;%深度
r=0.4; %频带带宽比例
alpha=0.2; %衰减系数
width=3.81;%宽度
c=1.54*10^5;%速度
d=(0:points-1).*depth/points;%每个点深度
fs=points*c/2/depth; %采样频率41.53MHz
fc=10000000;%探头中心频率
depth_pix=512;
width_pix=round(depth_pix/depth*width);
gap=fs/1618;
time=(0:points-1)/fs;
freq=0:gap:fs-fs/1618; %实际每个点频率
for line_no=1:lines
    a_line=img(line_no,:);
    number=66;
    %% 去除直流
    %高通滤波器设计
    wp=pi/2;ws=pi/7;
    highpass=wp-ws;
    N0=ceil(6.8*pi/highpass);
    N_highpass=N0+mod(N0+1,2);
    wc=(wp+ws)/2/pi; %中心频率
    highpass_filt=fir1(N_highpass-1,wc,'high',hanning(N_highpass)); %FIR butterworth highpass filter
    a_line_filtered=fftfilt(highpass_filt,a_line);
    figure(1); %原始及高通滤波
    if(line_no==number) %画一条看一看
        subplot(221);
        plot(time,a_line);title('原始射频时域');
        subplot(222);plot(time,a_line_filtered);title('滤波后时域');
        subplot 223;plot(freq,abs(fft(a_line)));title('原始射频频谱');
        subplot 224;plot(freq,abs(fft(a_line_filtered)));title('高通滤波后频谱');
    end
    
    %% 正交解调
    f_shift=0.11513*(fc*r)^2/(2*log(2))*alpha.*d;
    f_1=fc-f_shift; %动态频率
    Q=a_line_filtered.*cos(2*pi*f_1.*time);
    I=a_line_filtered.*sin(2*pi*f_1.*time);
    if(line_no==number) 
        figure(2); %正交解调
        subplot(221);
        plot(Q);%hold on;plot(a_line_filtered);
        title('解调后的Q');
        subplot(222);
        plot(I);%hold on;plot(a_line_filtered);
        title('解调后的I');
        subplot(223);
        plot(abs(fft(Q)));%hold on;plot(freq,abs(fft(a_line_filtered)));
        title('解调后的Q in 频域');
        subplot(224);
        plot(abs(fft(I)));%hold on;plot(freq,abs(fft(a_line_filtered)));
        title('解调后的I in 频域');
    end
    %% 高斯低通滤波
    lp_width=0.03;
    dz=depth/points;
    numz=2*round(2*lp_width/dz)+1;
    z2=dz*(-numz/2:(numz/2+1))';
    sigma=lp_width/4; %标准差
    LPF=((2*pi)^.5*sigma)^-1*exp(-.5*(z2/sigma).^2);
    Q_filtered=fftfilt(LPF,Q);
    I_filtered=fftfilt(LPF,I);
    lp_result=sqrt(Q_filtered.^2+I_filtered.^2);
    %画图看看
    if(line_no==number)
        figure(6);
        subplot 221
        plot(time,a_line_filtered);title('低通滤波前的时域');
        subplot 222
        plot(freq,abs(fft(a_line_filtered)));title('低通滤波前的频域');
        subplot 223
        plot(time,lp_result);title('低通滤波后的时域');
        subplot 224
        plot(freq,abs(fft(lp_result)));title('低通滤波后的频域');
    end
    
    %% 时间增益补偿
    beta=log(10)*alpha*fc/20;
    TGC_Matrix=1-exp(-beta.*d);
    I_TGC=I_filtered.*TGC_Matrix;
    Q_TGC=Q_filtered.*TGC_Matrix;
    envelop_IQ=sqrt(I_TGC.^2+Q_TGC.^2); %包络
    if(line_no==number)
        figure(3); %下采样
        subplot 221
        plot(abs(envelop_IQ));title('下采样前时域波形');
        subplot 222;
        plot(abs(fft(envelop_IQ)));title('下采样前频谱');
    end
    %% 下采样
    per_pix=fix(points/512);
    downsp=zeros(1,512);
    for i=1:512
        downsp(1,i)=envelop_IQ(1,3*(i-1)+1);
    end
    if(line_no==number)
        subplot 223;
        plot(abs(downsp));title('下采样后时域波形');
        subplot 224;
        plot(abs(fft(downsp)));title('下采样后频谱');
    end
    one_frame(line_no,:)=downsp;
end
% % %% 对数变换
% D=100;G=0;
% % for nlog=1:168
% %     for i=1:512
% %         q=D*log10(abs(one_frame(nlog,i))+1)+G;
% %         if q>255
% %             q=255;
% %         end
% %         logdata(nlog,i)=q;
% %     end
% % end
%% 对数变换new
D=60;G=0;
logdata=D*log10(one_frame+1)+G;
logdata=logdata';
figure;
imshow(one_frame',[]);
title('重建图像');
figure;
imshow(logdata,[]);
title('对数变换后的图像');
%% 插值显示
for i=1:size(logdata,1)
    for j=1:size(logdata,2)-1
        Interp_Out(i,2*j-1)=logdata(i,j);
        Interp_Out(i,2*j)=logdata(i,j)+0.5*(logdata(i,j+1)-logdata(i,j));
    end
end

figure 
b=(Interp_Out-min(Interp_Out(:)))./(max(Interp_Out(:))-min(Interp_Out(:)))*255;%%a为double型
b=b-50;
imshow(uint8(b));   
title('插值后的图像')
