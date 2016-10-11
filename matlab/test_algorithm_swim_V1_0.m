%Specification：游泳暂停，转向，划水次数检测            
%Reference：
%Author：xxf
%Time  ：2016-09-10
%Copyright：YF
%需要改善：滤波算法，谷峰检测算法，次数与圈数的不同步，2016-09-12
%缺   点：对泳池中途换泳姿计次不是很理想，圈数延迟问题（大概延迟半圈）
close all
clear
clc
% 01、C11 1085831：自由泳 7 7 7 7 7 7 :42
% 02、C12 1090302：蛙  泳 7 7 7 8  :29
% 03、C13 1090603：仰  泳 7 7 6 7  :27
% 04、C14 1091017：蝶  泳 8 9 7 8  :32
% 05、C21 1092628：蝶  泳 9 10 8 8  :35
% 06、C22 1093003：混合泳 8 8（蝶） 6 8（仰） 8 8（蛙） 7 7（自）:60
% 07、C23 1093911：自由泳 6 (4+3)（模拟暂停） 6 7  :26 圈数少1
% 08、C24 1094158：自由泳 7 7（模拟绕圈） 7 7  :28
% 09、C25 1094545：自由泳 7 7 7 7（右手）  :28
% 10、C31 1095154：自由泳（18圈）
% 11、C32 1095941：蛙  泳 6 7 7 7 7 7 7 7  :55
% 12、C33 1100354：仰  泳 6 6 6 6 6 7 7 8（前四圈可能各少记一下）:52
% 13、C34 1101036：混合泳 4 4 3 4（蝶->仰->蛙->自，中途换姿势）:15
% 14、C35 1101316：混合泳 4 4 3 4（蝶->仰->蛙->自，中途换姿势）:15
% 注：泳池长25米，表距离臂关节54cm，肘关节25cm，臂展182cm
selection = {'Mahony','data13'};%selection: 姿态算法(MG,Madgwick,Mahony)、游泳数据
% 01、C41：蛙  泳 9  9  :18 
% 02、C42：仰  泳 11 12-13 :23-24
% 03、C43：仰  泳 16 17 :33 p1奇葩
% 04、C44：蛙  泳 32 34 :66 p1
% 05、C45：蝶  泳 17 21 :39//
% 06、C46：自由泳 14 15 :29
% 07、C47：仰  泳 18 17 :35
% 08、C49：蛙  泳 11 12 :23 //
% 09、C410：混合泳 11 12 15 :38 圈数误差大
% 10、C411：蛙  泳 30 36 :66 p1
% 11、C51：自由泳 17 17 :34 p3
% 12、C52：仰  泳 16 18 :34
% 13、C53：蛙  泳 18 22 :40
% 14、C54：自由泳 16 20 :36
data = swim.select(selection);   %get data  
N = length(data);
Fs = 25; %sample frequence
lg = Fs * 0.6;
lg2 = Fs * 0.2;
lm = Fs / 5 * 8;
fifo4 = [0 -20 -40 1]';
fifo4t = [0 0 0 0]';
q = [1 0 0 0]';
% q = [0.7071 0.7071 0 0]';
fifot = zeros(3,1);
fifoi =zeros(3,1);
fifott = zeros(3,1);
fifoit =zeros(3,1);
fifott2 = zeros(3,1);
fifoit2 =zeros(3,1);
fifom = zeros(3,4);
fifoacosq1 = zeros(lg,1);
fifoacosq2 = zeros(lg,1);
fifoacosq1_2 = zeros(lg,1);
fifoacosq2_2 = zeros(lg,1);
fifonm1 = zeros(lm,4);%mx my mz nm
fifonm2 = zeros(lm,4);
KFg = [0 0]';
flag = 1;
step = 2;
stroke = 0;
lap = 0;
lapt = 0;
c5 = 5;
KFmt = [0 -5]';
KFmn = [0 0]';
Euler = zeros(N,3);
maxKFmn = 150;
maxdet = 0;
cstatic = 0;
%-----以下temp开头的变量用于绘图分析，可直接删除-----
tempdet = zeros(N,2);
tempfnm2 = zeros(N,5);
tempnm = zeros(N,1);
tempacosq = zeros(N,1);
tempacosq1 = zeros(N,1);
tempacosq2 = zeros(N,1);
tempacosq1_2 = zeros(N,1);
tempacosq2_2 = zeros(N,1);
tempacosq2KF = zeros(N,2);
tempacosq2KFt = zeros(N,2);
tempacosq2KFt2 = zeros(N,2);
temp_t_stroke_lap = zeros(N,4);
tempfifomKF = zeros(N,2);
tempdet16 = zeros(N,7);
tempdiffKFmn = zeros(N,2);
tempfifo4t = zeros(N,1);
%--------------------------------------------
A1 =[0.00272457781241851,-0.00153457263573283,-7.54004141204668e-05
    0,0.00192881941306130,-0.000390509025963946
    0,0,0.00190987873307018];
b1 = [69.2010347220676,-355.927726253649,-317.252668090672];

A2 = [0.00184660378338664,-0.000160362754257810,-0.000414085924444608
    0,0.00180115668441700,-0.000169701546512634
    0,0,0.00217324287849178];
b2 = [169.435864534901,-145.440551430268,-177.352724123477];
A3 = [0.00192497020496001,-2.12776513521997e-05,-4.78251520457670e-05
    0,0.00189741361107982,4.47275912233477e-05
    0,0,0.00196775945058362];
b3 = [103.373273840440,-387.867984202147,-135.189890733092];
fifokf = zeros(6*25,1);
fifostd = zeros(2*25,3);
tempfifostd = zeros(N,4);
tempm = zeros(N,4);
for i=1:1:N%1

%     q = [1 0 0 0]';
    
    gx1 = data(i,5) * (0.5 * (1.0 / Fs));
    gy1 = data(i,6) * (0.5 * (1.0 / Fs));
    gz1 = data(i,7) * (0.5 * (1.0 / Fs));
    A = [0   -gx1 -gy1 -gz1; 
         gx1  0    gz1 -gy1; 
         gy1 -gz1  0    gx1; 
         gz1  gy1 -gx1  0];
    q = q + A*q;
    q = q / sqrt(q(1)*q(1)+q(2)*q(2)+q(3)*q(3)+q(4)*q(4));
    
%     if data(i,1) == 50
%         q = [0.7071 0.7071 0 0]';
%     end
%     if q(1) < -0.1736 %超过120度预防越界
%         q = [1 0 0 0]';
%     end
    Euler(i,:) = swim.q2e(q);%欧拉角，roll，pitch，yaw
    acosq = 2*acos(q(1))*180/pi;
    fifokf = [fifokf(2:6*25,1);acosq];
    m0_4 = sum(fifokf(1:2*25,1).*fifokf(1+0.4*25:2.4*25,1))/50;
    m0_8 = sum(fifokf(1:2*25,1).*fifokf(1+0.8*25:2.8*25,1))/50;
    m1_6 = sum(fifokf(1:2*25,1).*fifokf(1+1.8*25:3.8*25,1))/50;
    tempm(i,:) = [data(i,1), m0_4, m0_8, m1_6];
    fifostd = [fifostd(2:50,:);[m0_4, m0_8, m1_6]];
    tempfifostd(i,:) = [data(i,1) sum(fifostd-repmat(mean(fifostd),50,1))/50];
%     acosq = acosq + gz1*180/pi;
    tempacosq(i,:) = acosq;
    
    fifoacosq1 = [fifoacosq1(2:lg,:);acosq];
    acosq1 = sum(fifoacosq1)/lg;
    tempacosq1(i,:) = acosq1;
    fifoacosq2 = [fifoacosq2(2:lg,:);acosq1];
    acosq2 = sum(fifoacosq2)/lg;
    tempacosq2(i,:) = acosq2;
    
    fifoacosq1_2 = [fifoacosq1_2(2:lg2,:);acosq];
    acosq1_2 = sum(fifoacosq1_2)/lg2;
    tempacosq1_2(i,:) = acosq1_2;
    fifoacosq2_2 = [fifoacosq2_2(2:lg2,:);acosq1_2];
    acosq2_2 = sum(fifoacosq2_2)/lg2;
    tempacosq2_2(i,:) = acosq2_2;
    
    fifott = [fifott(2:3,:);data(i,1)];
    fifoit = [fifoit(2:3,:);acosq];
    if((fifoit(3,1) >= fifoit(2,1) && fifoit(1,1) - fifoit(2,1) > 10e-7)||...
       (fifoit(3,1) <= fifoit(2,1) && fifoit(2,1) - fifoit(1,1) > 10e-7))
        tempacosq2KFt(i,1) = fifott(2,1);
        tempacosq2KFt(i,2) = fifoit(2,1);
    end
    fifott2 = [fifott2(2:3,:);data(i,1)];
    fifoit2 = [fifoit2(2:3,:);acosq2_2];
    if((fifoit2(3,1) >= fifoit2(2,1) && fifoit2(1,1) - fifoit2(2,1) > 10e-7)||...
       (fifoit2(3,1) <= fifoit2(2,1) && fifoit2(2,1) - fifoit2(1,1) > 10e-7))
        tempacosq2KFt2(i,1) = fifott2(2,1);
        tempacosq2KFt2(i,2) = fifoit2(2,1);
    end
    
    fifot = [fifot(2:3,:);data(i,1)];
    fifoi = [fifoi(2:3,:);acosq2];%滤波器选择
%     if fifoi(3,1) <= fifoi(2,1) && fifoi(2,1) - fifoi(1,1) > 10e-7
%        q = [1 0 0 0]';
%     end
    if((fifoi(3,1) >= fifoi(2,1) && fifoi(1,1) - fifoi(2,1) > 10e-7)||...
       (fifoi(3,1) <= fifoi(2,1) && fifoi(2,1) - fifoi(1,1) > 10e-7))
       tempacosq2KF(i,1) = fifot(2,1);
       tempacosq2KF(i,2) = fifoi(2,1);
       KFg = [KFg(2,1);fifoi(2,1)];
       det = KFg(2,1) - KFg(1,1);%------parameter
       tempdet(i,1) = data(i,1);
       tempdet(i,2) = abs(det);
       cstatic = cstatic + 1;
       if abs(det) > 15 || (abs(det)/abs(fifo4(4,1)) > 0.6 && abs(det) > 8 &&  abs(det)/maxdet > 0.1)||...
               (cstatic >= 4 && abs(det) > 8 && abs(det)/maxdet > 0.1)
           cstatic = 0; 
           fifo4t = [fifo4t(2:4,1); data(i,1)];
           fifo4 = [fifo4(2:4,1); det];
           tempdet16(i,1:2) = [data(i,1) abs(det)+20];
%        end
%       
           s1 = abs(abs(fifo4(1,1)) - abs(fifo4(2,1)));
           s2 = abs(abs(fifo4(3,1)) - abs(fifo4(4,1)));
           s3 = 2*s1/(abs(fifo4(1,1))+abs(fifo4(2,1)));%第一个三角形左腰与右腰的大小差
           s4 = 2*s2/(abs(fifo4(3,1))+abs(fifo4(4,1)));%第二个三角形左腰与右腰的大小差
           s5 = abs(abs(fifo4(1,1)) - abs(fifo4(3,1)));%第一个三角形与第二个三角形腰差
           tempdet16(i,3:7) = [s1 s2 s3 s4 s5];
           if step == 2
               if ((s3 < 0.18 && s4 < 0.18) || (s1 < 11 && s2 < 11)) && s5 < 48
                   if flag == 1
%                        disp(['t3: ',num2str(data(i,1)),' ','strokes: ',num2str(stroke),' ','lapt: ',num2str(lapt)]);
                       lap = lapt;
                       stroke = stroke + 3;
                       flag = 2;
                       tempstroke(i-1:i,:) = [fifo4t(3:4,:) fifo4(3:4,1)];
                   else
                       stroke = stroke + 1;
                       tempstroke(i-1:i,:) = [fifo4t(3:4,:) fifo4(3:4,:)];
                   end
                   if abs(det) > maxdet
                       maxdet = abs(det);
                   end
                   tempfifo4t(i,1) = fifo4t(4,1) - fifo4t(1,1);
                   step = 1;
               else
                   flag = 1;
%                    q = [1 0 0 0]';
%                    fifoacosq1 = zeros(lg,1);
               end
           else
               step = step + 1;
           end
       end
    end
    
    if c5 == 5
        c5 = 1;
%         tdata = AAA1*(data(i,8:10)-bbb1).';%加磁力计校准
        tdata = [data(i,8:10) 0];%不加磁力计校准
        nm = sqrt(tdata(:,1)*tdata(:,1) + tdata(:,2)*tdata(:,2) + tdata(:,3)*tdata(:,3));
        tempnm(i,:) = nm;
        tdata = [tdata(:,1:3) nm];
        
        fifonm1 = [fifonm1(2:lm,:);tdata];
        fnm1 = sum(fifonm1)/lm;
        fifonm2 = [fifonm2(2:lm,:);fnm1];
        fnm2 = sum(fifonm2)/lm;
        tempfnm2(i,:) = [data(i,1) fnm2];%滤波后数据
%         disp(['t: ',num2str(data(i,1)),' ','fmn2: ',num2str(fnm2)])
        fifot = [fifot(2:3,1);data(i,1)];
        fifom = [fifom(2:3,:);fnm2];
        ch_m = 2;
        if(((fifom(3,ch_m) >= fifom(2,ch_m)) && (fifom(1,ch_m) - fifom(2,ch_m) > 10e-7)) ||...
           ((fifom(3,ch_m) <= fifom(2,ch_m)) && (fifom(2,ch_m) - fifom(1,ch_m) > 10e-7)))
%             disp(['KF time: ', num2str(fifot(2,1))]);
            tempfifomKF(i,:) = [data(i,1) fifom(2,ch_m)];
            if fifot(2,1) - KFmt(2,1) >= 5 && stroke > 0%峰谷时间差大于5s，次数大于0
                KFmt = [KFmt(2,1); fifot(2,1)];
                KFmn = [KFmn(2,1); fifom(2,ch_m)];%移位存储新的谷峰点
                diffKFmn = abs(KFmn(2,1) - KFmn(1,1));
                tempdiffKFmn(i,:) = [data(i,1) diffKFmn];
                if diffKFmn > maxKFmn %初始值150
                    maxKFmn = diffKFmn;
                end
                dd = diffKFmn/maxKFmn;

                if (diffKFmn > 170 && dd > 0.55) || (dd > 0.9 && lapt == 0)||diffKFmn > 200%主要判别依据：磁力计波动大于180uT
                    lapt = lapt + 1;  
%                     disp(['t1: ',num2str(data(i,1)),' ','strokes: ',num2str(stroke),' ','lapt: ',num2str(lapt)]);
                elseif lapt == 0
                    lapt = lapt + 1;
%                     disp(['t2: ',num2str(data(i,1)),' ','strokes: ',num2str(stroke),' ','lapt: ',num2str(lapt)]);
                end
                disp(['t: ',num2str(KFmt(2,1)),' ','strokes: ',num2str(stroke),' ','lapt: ',num2str(lapt)]);
            end
        end
    else
        c5 = c5 + 1;
    end
%     disp(['t: ',num2str(data(i,1)),' ','strokes: ',num2str(stroke),' ','laps: ',num2str(lap)]);
    temp_t_stroke_lap(i,:) = [data(i,1), stroke, lap, lapt];
end
disp(['all strokes: ',num2str(stroke),'   ','all laps: ',num2str(lapt)]);
%explot(data,extitle,excolor,exlabel,exlegend,exline)
%explot2(data,row,extitle,excolor,exlabel,exlegend,exline)
tempfifo4t(tempfifo4t(:,1) == 0,:) = [];
% av = sum(tempfifo4t)/length(tempfifo4t)*2/3
% ma = max(tempfifo4t)*2/3
%--------------原始数据-------------------
figure(1)%欧拉角，三轴加速度、陀螺仪、磁力计
swim.plotoriginal(Euler,data,selection(2),'data of original-')%,staticdata)
%---------------次数-----------------------
figure(2)%acosq
tempacosq2KF(tempacosq2KF(:,2) == 0,:) = [];
tempacosq2KFt(tempacosq2KFt(:,2) == 0,:) = [];
subplot(5,1,1)
plot(data(:,1),tempacosq)
hold on
plot(tempacosq2KFt(:,1),tempacosq2KFt(:,2),'r*')
title('original-acosq')
xlabel('t(s)');ylabel('degree(°)');
legend('acosq','KF of acosq')
hold on 
plot(data(:,1),360*ones(size(data(:,1))),'k','linewidth',1)
subplot(5,1,2)
plot(data(:,1),tempacosq2_2)
hold on
plot(tempacosq2KFt2(:,1),tempacosq2KFt2(:,2),'r*')
title('acosq-0.2')
xlabel('t(s)');ylabel('degree(°)');
legend('acosq','KF of acosq')
subplot(5,1,3)
plot(data(:,1),tempacosq2)
hold on
plot(tempacosq2KF(:,1),tempacosq2KF(:,2),'r*')
title('acosq-long')
xlabel('t(s)');ylabel('degree(°)');
legend('acosq','KF of acosq')
subplot(5,1,4)%det
tempdet(tempdet(:,2) == 0,:) = [];
tempdet16(tempdet16(:,3) == 0,:) = [];
stem(tempdet(:,1),abs(tempdet(:,2)))
hold on
stem(tempdet16(:,1),tempdet16(:,2),'g');
hold on
stem(tempstroke(:,1),abs(tempstroke(:,2)),'r')
hold on
plot(tempdet(:,1),15*ones(size(tempdet(:,1))),'k')
legend('original','valid','stroke','threshold')
title('difference of acosq')
xlabel('t(s)');ylabel('degree(°)');
subplot(5,1,5)
swim.explot(temp_t_stroke_lap,'t-stroke-lap',{'r','g','b'},{'t(s)','times'},{'strokes','laps','lapt'},{'-','-','-'})
grid on
% ----------------圈数-----------------------
% figure(3)%mag
% tempfnm2(tempfnm2(:,2) == 0,:) = [];
% tempnm(tempnm(:,1) == 0,:) = [];
% suptitle('未加校准磁力计')
% subplot(3,1,1)
% swim.explot2([[data(:,1),data(:,9)];[tempfnm2(:,1),tempfnm2(:,3)];tempfifomKF;...
%     [temp_t_stroke_lap(:,1),200*temp_t_stroke_lap(:,4)];[temp_t_stroke_lap(:,1),200*temp_t_stroke_lap(:,3)]],...
%     [length(data(:,1)),length(tempfnm2(:,1)),length(tempfifomKF),length(temp_t_stroke_lap),length(temp_t_stroke_lap)],...
%     'my',{'g','k','r','m','b'},{'t(s)','uT'},{'before filter','after filter','KF','200*lapt-orig','200*laps-real'},{'-','-','*','-','-'});
% subplot(3,1,2)
% stem(tempdiffKFmn(:,1),tempdiffKFmn(:,2))
% hold on
% plot(0:data(end,1),170*ones(size(0:data(end,1))))
% title('difference of mag')
% xlabel('t(s)');ylabel('uT');
% legend('diff mag','threshold')
% subplot(3,1,3)
% swim.explot(temp_t_stroke_lap,'t-stroke-lap',{'r','b','m'},{'t(s)','times'},{'strokes','laps','lapt'},{'-','-','-'})
% grid on
% figure(4)
% plot(tempm(:,1),tempm(:,2))
% hold on
% plot(tempm(:,1),tempm(:,3))
% hold on
% plot(tempm(:,1),tempm(:,4))
% legend('m1-0','m2-0','m3-0')
% figure(5)
% 
% plot(tempfifostd(:,1),tempfifostd(:,2))
% hold on
% plot(tempfifostd(:,1),tempfifostd(:,3),'linewidth',1)
% hold on
% plot(tempfifostd(:,1),tempfifostd(:,4))
% legend('m1-0','m2-0','m3-0')


