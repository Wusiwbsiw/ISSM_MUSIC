clc;
clear ;
close all;

c = 3e8;

M = 8;%阵元数
K = 2;%时域信号段数
N = 256;%每段信号时域快照数
J = 256;%每段时域信号DFT点数
fs = 400e6;
fc = 100e6;
lambda=c/fc;
d=lambda/2;
B = 40e6;

SNR = 50;

P = 2;%信源数
target_angle = [19 ; 30] /180*pi;%信源角度

t = 0:1/fs:K*N/fs-1/fs;

Signal_transmit = [ chirp(t,fc-B/2,K*N/fs-1/fs,fc+B/2) ];
 Signal_transmit=awgn(Signal_transmit,SNR,'measured');              %在信号中添加高斯噪声

f = linspace(-fs/2,fs/2,length(Signal_transmit));
figure;plot(f,abs(fftshift(fft(Signal_transmit))));
figure;plot(t,Signal_transmit);

% 将信号分段
Segmented = zeros(K,N);
for i = 1:K
    Segmented(i,:) = Signal_transmit(1+(i-1)*N:i*N);
end
% 分段后DFT
DFTed = fftshift(fft(Segmented,J,2),2);
DFTed1 = DFTed;DFTed2 = DFTed;
figure;imagesc(f,1:K,abs(DFTed));
% 得到每个频点的协方差矩阵
RX = zeros(M,M,J);
f_vector =  linspace(fc-B/2,fc+B/2,J);%linspace(-fs/2,fs/2,J);%linspace(fc-B/2,fs/2,J);%
sump_music=zeros(1,361);       %角度功率谱初始化

figure;hold on;
for j = 1:J
    a = exp(-1j*2*pi*f_vector(j)*d*sin(target_angle)/c*(0:M-1));
    RX(:,:,j) = ([DFTed1(:,j) DFTed2(:,j)]*a)'*([DFTed1(:,j) DFTed2(:,j)]*a);
    %RX(:,:,j) = (ones(1,2)*a)'*(ones(1,2)*a);
    % 对每个频点的协方差矩阵进行分解估计方向
    [Ev,D] = eig(RX(:,:,j));     % 特征值分解 D：特征值的对角矩阵 Ev：右特征列向量组成的矩阵
    EVA = diag(D)';            % 将特征值提取为1行
    [EVA,I] = sort(EVA);       % 对特征值排序,从小到大。其中I为索引向量
    EV = fliplr(Ev(:,I));      % 按照索引I对顺序特征矢量排序得到Ev,再fliplr水平颠倒列向量得到特征值从大到小分布的特征列向量组成的矩阵EV
    En = EV(:,P+1:end);          % 取特征向量矩阵的第X+1到M列特征向量组成噪声子空间En
    for i = 1:361
        angle(i) = (i-181)/2;           % 映射到-90度到90度
        theta_m = angle(i)*pi/180;
        a_theta = exp(-1j*2*pi*f_vector(j)/c*d*(0:M-1)*sin(theta_m));    %导向矢量M*1
        p_music(i) = abs(1/(a_theta*En*En'*a_theta'));          %MUSIC算法功率谱
    end
    sump_music=sump_music+p_music;                              %累加各频点功率谱
    plot(angle,p_music);
end



%%ISM_MUSIC
p_music_max = max(sump_music);
sump_music = 10*log10(sump_music/p_music_max);
figure;
plot(angle,sump_music,'b-');
%title('ISM——MUSIC空间谱');
xlabel('入射角/度');
ylabel('空间谱/dB');