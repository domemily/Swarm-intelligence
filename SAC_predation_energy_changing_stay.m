clc
clear all
rand('seed',1);
%% 本程序主要修改内容如下：
% 1. 进食时间改为固定值
% 2. 噪音加在最终速度项
% 3. 
%% 初始化参数
n_1 = 100;    %被捕食者的数目
n_2 = 5;   %捕食者的数目
timestep=1000;
delta=0.05;
delay = 50; %进食时间

v_co = 8; 
v_to = 6; 
r_0_1 = 0;  %prey 的速度协同半径
r_e_1 = 1;  %prey 的位置协同半径
r_0_2 = 0;  %pred 的速度协同半径 
r_e_2 = 1;  %pred 的位置协同半径
r_s = 10;  %捕食交互半径
r_c = 1; %被捕半径
d0=1;   %个体间最小间距
L = 40; %仿真空间大小
noise_pred = 0.1; %捕食者速度噪声
noise_prey = 0.1; %猎物速度噪声
noise_active = 1; %捕食者再次捕食的方向随机量

alfa = 1; %速度匹配强度
beta=5;   %位置协同强度 
gamma=1.2;  %追捕强度

k_1 = 1;  %没有天敌时的能量递增系数
k_2 = 1;  % 同伴牺牲时的群体收益系数
k_3 = -0.5; %逃逸时猎物的能量消耗系数
k_4 = -0.1; % 追捕时捕食者的能量消耗系数

unit_e = 0.01; %能量递增/递减单元
in_e = 0.8; %初始能量值，最大能量值为100；
max_e = 1; %最大能量值
min_e = 0.1;
z_pred=zeros(n_2,5,timestep); %状态量声明
z_prey=zeros(n_1,5,timestep); 
shap_1=zeros(n_1,2,timestep);
shap_2=zeros(n_1,2,timestep);
shap_3=zeros(n_1,2,timestep);
shap_4=zeros(n_1,2,timestep);
shap_1_d=zeros(n_2,2,timestep);
shap_2_d=zeros(n_2,2,timestep);
shap_3_d=zeros(n_2,2,timestep);
shap_4_d=zeros(n_2,2,timestep);
%% 给定随机的初始位置
a=2*L*rand(n_1,4)-L;b=2*L*rand(n_2,4)-L;c=in_e*ones(n_1,1);d=in_e*ones(n_2,1);
z_pred(:,:,1)=[b,d];
z_prey(:,:,1)=[a,c];
%% 初始化结构体
pos_datax_22=zeros(n_2,n_2);
pos_datay_22=zeros(n_2,n_2);
spe_datax_22=zeros(n_2,n_2);
spe_datay_22=zeros(n_2,n_2);
pos_datax_11=zeros(n_1,n_1);
pos_datay_11=zeros(n_1,n_1);
spe_datax_11=zeros(n_1,n_1);
spe_datay_11=zeros(n_1,n_1);
datax_11=zeros(n_1,n_1);
datay_11=zeros(n_1,n_1);
datax_pred=zeros(n_2,n_1);
datay_pred=zeros(n_2,n_1);
chase_21=zeros(n_2,n_1);
chase_12=zeros(n_1,n_2);
f_ct_x=zeros(n_2,1);
f_ct_y=zeros(n_2,1);
prey_worse=zeros(n_1,n_2);
prey_normal=zeros(n_1,n_2);
stop_pred=zeros(n_2,timestep);

gix_2=zeros(n_2,1);
giy_2=zeros(n_2,1);
pix_2=zeros(n_2,1);
piy_2=zeros(n_2,1);
gix_1=zeros(n_1,1);
giy_1=zeros(n_1,1);
pix_1=zeros(n_1,1);
piy_1=zeros(n_1,1);

%% 初始化位置和速度
gix_2(:)=z_pred(:,1,1);
giy_2(:)=z_pred(:,2,1);
pix_2(:)=z_pred(:,3,1);
piy_2(:)=z_pred(:,4,1);
gix_1(:)=z_prey(:,1,1);
giy_1(:)=z_prey(:,2,1);
pix_1(:)=z_prey(:,3,1);
piy_1(:)=z_prey(:,4,1);
 %% 初始时刻捕食者和猎物的能量总和：
 total_energy_pred(1)=n_2*in_e;
 total_energy_prey(1)=n_1*in_e;
 nature_energy(1)=total_energy_pred(1)+total_energy_prey(1);
 num_live_prey(1)=100;
%% 程序进入主循环
for jj=2:1:timestep
    %% 捕食者更新程序
   for ii=1:n_2
    for k=1:n_2
          % 邻居定义为平面坐标距离小于d的个体：
          d_22=sqrt((z_pred(k,1,jj-1)-z_pred(ii,1,jj-1))^2+...
              (z_pred(k,2,jj-1)-z_pred(ii,2,jj-1))^2);
          if d_22<=r_e_2  %计算位置协同项         
            if d_22==0;
            pos_datax_22(ii,k)=0;
            pos_datay_22(ii,k)=0; 
            else
            valx=2*((1/d_22)-(d0/(d_22^3)))*(z_pred(k,1,jj-1)-z_pred(ii,1,jj-1))/d_22; %近距离排斥远距离吸引
            valy=2*((1/d_22)-(d0/(d_22^3)))*(z_pred(k,2,jj-1)-z_pred(ii,2,jj-1))/d_22;
            pos_datax_22(ii,k)=valx;
            pos_datay_22(ii,k)=valy;
            end
          else
            pos_datax_22(ii,k)=0;
            pos_datay_22(ii,k)=0;
          end
          if d_22<=r_0_2  %计算速度协同项
            vvalx=z_pred(k,3,jj-1);
            vvaly=z_pred(k,4,jj-1);
            spe_datax_22(ii,k)=vvalx;
            spe_datay_22(ii,k)=vvaly;
          else
            spe_datax_22(ii,k)=0;
            spe_datay_22(ii,k)=0;
          end
    
    end
    pos_x_22=sum(pos_datax_22,2);
    pos_y_22=sum(pos_datay_22,2);
    spe_x_22=sum(spe_datax_22,2);
    spe_y_22=sum(spe_datay_22,2);
    %% 速度加噪声
   %  [x*cosA-y*sinA  x*sinA+y*cosA] 
    %% 求最小猎物
    
    %% 求捕食者与猎物的交互关系
    
        val_a=z_prey(:,1,jj-1);val_b=z_prey(:,2,jj-1);val_l=L*ones(n_1,1);
        shap_1(:,:,jj-1)=[val_a+val_l,val_b];shap_2(:,:,jj-1)=[val_a-val_l,val_b];
        shap_3(:,:,jj-1)=[val_a,val_b+val_l];shap_4(:,:,jj-1)=[val_a,val_b-val_l];
    for kk=1:n_1
        if z_prey(kk,5,jj-1)==-1000  %如果猎物已经死亡
           chase_21(ii,kk)=1000;  % 不再影响后续计算
           E_ind(ii,kk)=0; datax_pred(ii,kk) = 0; datay_pred(ii,kk) = 0; 
        else   % 计算各猎物在四个映射区里距离捕食者的距离
            d_21=sqrt((z_pred(ii,1,jj-1)-z_prey(kk,1,jj-1))^2+...
              (z_pred(ii,2,jj-1)-z_prey(kk,2,jj-1))^2);
            d_shap_1=sqrt((z_pred(ii,1,jj-1)-shap_1(kk,1,jj-1))^2+...
              (z_pred(ii,2,jj-1)-shap_1(kk,2,jj-1))^2);
            d_shap_2=sqrt((z_pred(ii,1,jj-1)-shap_2(kk,1,jj-1))^2+...
              (z_pred(ii,2,jj-1)-shap_2(kk,2,jj-1))^2);
            d_shap_3=sqrt((z_pred(ii,1,jj-1)-shap_3(kk,1,jj-1))^2+...
              (z_pred(ii,2,jj-1)-shap_3(kk,2,jj-1))^2);
            d_shap_4=sqrt((z_pred(ii,1,jj-1)-shap_4(kk,1,jj-1))^2+...
              (z_pred(ii,2,jj-1)-shap_4(kk,2,jj-1))^2);
          if d_21<=r_s && d_21>=r_c         %如果猎物在视野范围内，且在进食范围外
              val = d_21;
              chase_21(ii,kk)=val;
              valx = 1/(chase_21(ii,kk))^10*...
                     (z_pred(ii,1,jj-1)-z_prey(kk,1,jj-1));
              valy = 1/(chase_21(ii,kk))^10*...
                     (z_pred(ii,2,jj-1)-z_prey(kk,2,jj-1));
              datax_pred(ii,kk)=valx;
              datay_pred(ii,kk)=valy;  
              E_ind(ii,kk) = 0;
          elseif d_shap_1<=r_s && d_shap_1>=r_c         %如果猎物在视野范围内，且在进食范围外
              val = d_shap_1;
              chase_21(ii,kk)=val;
              valx = 1/(chase_21(ii,kk))^10*...
                     (z_pred(ii,1,jj-1)-shap_1(kk,1,jj-1));
              valy = 1/(chase_21(ii,kk))^10*...
                     (z_pred(ii,2,jj-1)-shap_1(kk,2,jj-1));
              datax_pred(ii,kk)=valx;
              datay_pred(ii,kk)=valy;  
              E_ind(ii,kk) = 0;
           elseif d_shap_2<=r_s && d_shap_2>=r_c         %如果猎物在视野范围内，且在进食范围外
              val = d_shap_2;
              chase_21(ii,kk)=val;
              valx = 1/(chase_21(ii,kk))^10*...
                     (z_pred(ii,1,jj-1)-shap_2(kk,1,jj-1));
              valy = 1/(chase_21(ii,kk))^10*...
                     (z_pred(ii,2,jj-1)-shap_2(kk,2,jj-1));
              datax_pred(ii,kk)=valx;
              datay_pred(ii,kk)=valy;  
              E_ind(ii,kk) = 0;
           elseif d_shap_3<=r_s && d_shap_3>=r_c         %如果猎物在视野范围内，且在进食范围外
              val = d_shap_3;
              chase_21(ii,kk)=val;
              valx = 1/(chase_21(ii,kk))^10*...
                     (z_pred(ii,1,jj-1)-shap_3(kk,1,jj-1));
              valy = 1/(chase_21(ii,kk))^10*...
                     (z_pred(ii,2,jj-1)-shap_3(kk,2,jj-1));
              datax_pred(ii,kk)=valx;
              datay_pred(ii,kk)=valy;  
              E_ind(ii,kk) = 0;
           elseif d_shap_4<=r_s && d_shap_4>=r_c         %如果猎物在视野范围内，且在进食范围外
              val = d_shap_4;
              chase_21(ii,kk)=val;
              valx = 1/(chase_21(ii,kk))^10*...
                     (z_pred(ii,1,jj-1)-shap_4(kk,1,jj-1));
              valy = 1/(chase_21(ii,kk))^10*...
                     (z_pred(ii,2,jj-1)-shap_4(kk,2,jj-1));
              datax_pred(ii,kk)=valx;
              datay_pred(ii,kk)=valy;  
              E_ind(ii,kk) = 0;
           elseif d_21<r_c
              E_ind(ii,kk) = z_prey(kk,5,jj-1); %如果进食区域有猎物，猎物能量输送给捕食者
              z_prey(kk,5,jj)=-1000; 
              chase_21(ii,kk)=1000;
              valx = 1/(chase_21(ii,kk))^10*...
                     (z_pred(ii,1,jj-1)-z_prey(kk,1,jj-1));
              valy = 1/(chase_21(ii,kk))^10*...
                     (z_pred(ii,2,jj-1)-z_prey(kk,2,jj-1));
              datax_pred(ii,kk)=valx;
              datay_pred(ii,kk)=valy; 
           elseif d_shap_1<r_c
              E_ind(ii,kk) = z_prey(kk,5,jj-1); %如果进食区域有猎物，猎物能量输送给捕食者
              z_prey(kk,5,jj)=-1000; 
              chase_21(ii,kk)=1000;
              valx = 1/(chase_21(ii,kk))^10*...
                     (z_pred(ii,1,jj-1)-shap_1(kk,1,jj-1));
              valy = 1/(chase_21(ii,kk))^10*...
                     (z_pred(ii,2,jj-1)-shap_1(kk,2,jj-1));
              datax_pred(ii,kk)=valx;
              datay_pred(ii,kk)=valy;
           elseif d_shap_2<r_c
              E_ind(ii,kk) = z_prey(kk,5,jj-1); %如果进食区域有猎物，猎物能量输送给捕食者
              z_prey(kk,5,jj)=-1000; 
              chase_21(ii,kk)=1000;
              valx = 1/(chase_21(ii,kk))^10*...
                     (z_pred(ii,1,jj-1)-shap_2(kk,1,jj-1));
              valy = 1/(chase_21(ii,kk))^10*...
                     (z_pred(ii,2,jj-1)-shap_2(kk,2,jj-1));
              datax_pred(ii,kk)=valx;
              datay_pred(ii,kk)=valy;
           elseif d_shap_3<r_c
              E_ind(ii,kk) = z_prey(kk,5,jj-1); %如果进食区域有猎物，猎物能量输送给捕食者
              z_prey(kk,5,jj)=-1000; 
              chase_21(ii,kk)=1000;
              valx = 1/(chase_21(ii,kk))^10*...
                     (z_pred(ii,1,jj-1)-shap_3(kk,1,jj-1));
              valy = 1/(chase_21(ii,kk))^10*...
                     (z_pred(ii,2,jj-1)-shap_3(kk,2,jj-1));
              datax_pred(ii,kk)=valx;
              datay_pred(ii,kk)=valy;
          elseif d_shap_4<r_c
              E_ind(ii,kk) = z_prey(kk,5,jj-1); %如果进食区域有猎物，猎物能量输送给捕食者
              z_prey(kk,5,jj)=-1000; 
              chase_21(ii,kk)=1000;
              valx = 1/(chase_21(ii,kk))^10*...
                     (z_pred(ii,1,jj-1)-shap_4(kk,1,jj-1));
              valy = 1/(chase_21(ii,kk))^10*...
                     (z_pred(ii,2,jj-1)-shap_4(kk,2,jj-1));
              datax_pred(ii,kk)=valx;
              datay_pred(ii,kk)=valy;
          else
              chase_21(ii,kk)=1000;
              valx = 1/(chase_21(ii,kk))^10;
              valy = 1/(chase_21(ii,kk))^10;
              datax_pred(ii,kk)=valx;
              datay_pred(ii,kk)=valy;
              E_ind(ii,kk) = 0;
          end
        end
    end
          ct_x(ii)=sum(datax_pred(ii,:));
          ct_y(ii)=sum(datay_pred(ii,:));
          f_ct_x(ii)=ct_x(ii)/sqrt((ct_x(ii))^2+(ct_y(ii))^2);  % 对附近猎物按距离加权求和
          f_ct_y(ii)=ct_y(ii)/sqrt((ct_x(ii))^2+(ct_y(ii))^2);
          
%% 捕食者速度位置更新方程
          uix_2(ii)=beta*pos_x_22(ii)+alfa*spe_x_22(ii)-gamma*f_ct_x(ii);
          uiy_2(ii)=beta*pos_y_22(ii)+alfa*spe_y_22(ii)-gamma*f_ct_y(ii);
          v_c(ii) = v_co*sigmf(z_pred(ii,5,jj-1),[10 0.5]); %速度更新方程
          
          pix_2(ii)=v_c(ii)*uix_2(ii)/(sqrt(uix_2(ii)^2+uiy_2(ii)^2));
          piy_2(ii)=v_c(ii)*uiy_2(ii)/(sqrt(uix_2(ii)^2+uiy_2(ii)^2));
          % 给速度增加噪音
          xita_pred = 2*noise_pred*pi*rand()-noise_pred*pi;
          pix_2(ii)=pix_2(ii)*cos(xita_pred)-piy_2(ii)*sin(xita_pred);
          piy_2(ii)=piy_2(ii)*cos(xita_pred)+pix_2(ii)*sin(xita_pred);

          gix_2(ii)=gix_2(ii)+pix_2(ii)*delta;
          giy_2(ii)=giy_2(ii)+piy_2(ii)*delta;
%% 能量更新程序
dot(ii) = z_pred(ii,3,jj-1)*pix_2(ii)+z_pred(ii,4,jj-1)*piy_2(ii);
module(ii) = sqrt((z_pred(ii,3,jj-1))^2+(z_pred(ii,4,jj-1))^2)*...
    sqrt((pix_2(ii))^2+(piy_2(ii))^2);
yuxuan(ii)=dot(ii)/module(ii);
theta(ii) = acos(yuxuan(ii));


E_add(ii) = sum(E_ind(ii,:));%将所有掠夺能量加和
z_pred(ii,5,jj) = z_pred(ii,5,jj-1) + E_add(ii) +...
    k_4*((v_c(ii)/v_co)^2+(theta(ii)/pi)^2)*unit_e; %能量更新方式
%% 捕食者能量存在上限，不需要下限 
%捕食者能量存在上限，但是不存在下限，能量取负值时，速度仍然为0，对程序不影响。 
if z_pred(ii,5,jj)>=max_e
    z_pred(ii,5,jj)=max_e;
end
%% 自由空间活动    
%     z_pred(ii,1,jj)=gix_2(ii);
%     z_pred(ii,2,jj)=giy_2(ii);
%     z_pred(ii,3,jj)=pix_2(ii);
%     z_pred(ii,4,jj)=piy_2(ii);
%% 捕食者的活动范围设置为周期环境
% 设置幅度    
if abs(gix_2(ii))>L/2
    gix_2(ii)=-sign(gix_2(ii))*L+gix_2(ii);
end
if abs(giy_2(ii))>L/2
    giy_2(ii)=-sign(giy_2(ii))*L+giy_2(ii);
end 

%% 捕食者状态量更新方程
    z_pred(ii,1,jj)=gix_2(ii);
    z_pred(ii,2,jj)=giy_2(ii);
    z_pred(ii,3,jj)=pix_2(ii);
    z_pred(ii,4,jj)=piy_2(ii);
%% 设置进食延迟
if jj>3
if z_pred(ii,5,jj-1)>z_pred(ii,5,jj-2) 
    stop_pred(ii,jj)=jj-1;
elseif z_pred(ii,5,jj-1)==z_pred(ii,5,jj-2) 
    stop_pred(ii,jj)=stop_pred(ii,jj-1);
else
    stop_pred(ii,jj)=jj-1-delay;
end
if jj<stop_pred(ii,jj)+delay
    z_pred(ii,:,jj)=z_pred(ii,:,jj-1);
    gix_2(ii) = z_pred(ii,1,jj); % gix和giy是中间变量，直接决定z，
    giy_2(ii) = z_pred(ii,2,jj);  % 容易只对z处理而忘了实际参与运算的是gix，giy
    
    % 下次捕食的方向随机产生，从而与独立生存的个体相匹配
    xita_active= 2*noise_active*pi*rand()-noise_active*pi;
    pix_2(ii)=pix_2(ii)*cos(xita_active)-piy_2(ii)*sin(xita_active);
    piy_2(ii)=piy_2(ii)*cos(xita_active)+pix_2(ii)*sin(xita_active);
    z_pred(ii,3,jj)=pix_2(ii);
    z_pred(ii,4,jj)=piy_2(ii);
    
    % 下次捕食时的位置随机产生，从而模拟数目不变，周期性的捕食行为。
    
   
end
end

   end
  
%% 逃逸者更新程序---------------------------------------------------------
    for ii=1:n_1
        if z_prey(ii,5,jj-1)==-1000          % 如果猎物被捕杀死亡，则不再参与运算
               z_prey(ii,1,jj)=NaN;      % 猎物将从仿真环境中移除
               z_prey(ii,2,jj)=NaN;
               z_prey(ii,5,jj)=-1000;        % 猎物能量被强制设为-10000
        else
    for k=1:n_1
          % 邻居定义为平面坐标距离小于r_s的个体：
          d_11=sqrt((z_prey(k,1,jj-1)-z_prey(ii,1,jj-1))^2+...
              (z_prey(k,2,jj-1)-z_prey(ii,2,jj-1))^2);
          if d_11<=r_e_1
            %计算位置协同项
            if d_11==0
                pos_datax_11(ii,k)=0;
                pos_datay_11(ii,k)=0;
            else
            valx=2*((1/d_11)-(d0/(d_11^3)))*(z_prey(k,1,jj-1)-z_prey(ii,1,jj-1))/d_11; %近距离排斥远距离吸引
            valy=2*((1/d_11)-(d0/(d_11^3)))*(z_prey(k,2,jj-1)-z_prey(ii,2,jj-1))/d_11;
            pos_datax_11(ii,k)=valx;
            pos_datay_11(ii,k)=valy;
            end
          else
            pos_datax_11(ii,k)=0;
            pos_datay_11(ii,k)=0;
          end
          if d_11<=r_0_1 %计算速度协同项
            vvalx=z_prey(k,3,jj-1);
            vvaly=z_prey(k,4,jj-1);
            spe_datax_11(ii,k)=vvalx;
            spe_datay_11(ii,k)=vvaly;
          else
            spe_datax_11(ii,k)=0;
            spe_datay_11(ii,k)=0;
          end
    end
    pos_x_11=sum(pos_datax_11,2);
    pos_y_11=sum(pos_datay_11,2);
    spe_x_11=sum(spe_datax_11,2);
    spe_y_11=sum(spe_datay_11,2);
    
 %% 求最近捕食者
   val_c=z_pred(:,1,jj-1);val_d=z_pred(:,2,jj-1);val_l_d=L*ones(n_2,1);
    shap_1_d(:,:,jj-1)=[val_c+val_l_d,val_d];shap_2_d(:,:,jj-1)=[val_c-val_l_d,val_d];
    shap_3_d(:,:,jj-1)=[val_c,val_d+val_l_d];shap_4_d(:,:,jj-1)=[val_c,val_d-val_l_d];  
    for kk=1:n_2
            d_12=sqrt((z_prey(ii,1,jj-1)-z_pred(kk,1,jj-1))^2+...
              (z_prey(ii,2,jj-1)-z_pred(kk,2,jj-1))^2);
            d_shap_1_d=sqrt((z_prey(ii,1,jj-1)-shap_1_d(kk,1,jj-1))^2+...
              (z_prey(ii,2,jj-1)-shap_1_d(kk,2,jj-1))^2);
            d_shap_2_d=sqrt((z_prey(ii,1,jj-1)-shap_2_d(kk,1,jj-1))^2+...
              (z_prey(ii,2,jj-1)-shap_2_d(kk,2,jj-1))^2);
            d_shap_3_d=sqrt((z_prey(ii,1,jj-1)-shap_3_d(kk,1,jj-1))^2+...
              (z_prey(ii,2,jj-1)-shap_3_d(kk,2,jj-1))^2);
            d_shap_4_d=sqrt((z_prey(ii,1,jj-1)-shap_4_d(kk,1,jj-1))^2+...
              (z_prey(ii,2,jj-1)-shap_4_d(kk,2,jj-1))^2);
          if d_12<=r_s && d_12>=r_c  %视野内有捕食者，但自己没被捕杀
              val = d_12;
              chase_12(ii,kk)=val;
              valx = 1/(chase_12(ii,kk)-r_c)^2*...
                     (z_prey(ii,1,jj-1)-z_pred(kk,1,jj-1));
              valy = 1/(chase_12(ii,kk)-r_c)^2*...
                     (z_prey(ii,2,jj-1)-z_pred(kk,2,jj-1));
              datax(ii,kk)=valx;
              datay(ii,kk)=valy;  
              prey_normal(ii,kk)=0;
              if jj>3
              if z_pred(kk,5,jj-1)>=z_pred(kk,5,jj-2)  % 这里注意对捕食行为的判断
                  prey_worse(ii,kk)=0;        %是建立在能量不减少，或者能量增加的；
                                              %如果捕食者达到能量下限的情况不变，则会被误认为进食。
              else
                  prey_worse(ii,kk)=1;

              end
              end
          elseif d_shap_1_d<=r_s && d_shap_1_d>=r_c  %视野内有捕食者，但自己没被捕杀
              val = d_shap_1_d;
              chase_12(ii,kk)=val;
              valx = 1/(chase_12(ii,kk)-r_c)^2*...
                     (z_prey(ii,1,jj-1)-shap_1_d(kk,1,jj-1));
              valy = 1/(chase_12(ii,kk)-r_c)^2*...
                     (z_prey(ii,2,jj-1)-shap_1_d(kk,2,jj-1));
              datax(ii,kk)=valx;
              datay(ii,kk)=valy;  
              prey_normal(ii,kk)=0;
              if jj>3
              if z_pred(kk,5,jj-1)>=z_pred(kk,5,jj-2)  % 这里注意对捕食行为的判断
                  prey_worse(ii,kk)=0;        %是建立在能量不减少，或者能量增加的；
                                              %如果捕食者达到能量下限的情况不变，则会被误认为进食。
              else
                  prey_worse(ii,kk)=1;

              end
              end
          elseif d_shap_2_d<=r_s && d_shap_2_d>=r_c  %视野内有捕食者，但自己没被捕杀
              val = d_shap_2_d;
              chase_12(ii,kk)=val;
              valx = 1/(chase_12(ii,kk)-r_c)^2*...
                     (z_prey(ii,1,jj-1)-shap_2_d(kk,1,jj-1));
              valy = 1/(chase_12(ii,kk)-r_c)^2*...
                     (z_prey(ii,2,jj-1)-shap_2_d(kk,2,jj-1));
              datax(ii,kk)=valx;
              datay(ii,kk)=valy;  
              prey_normal(ii,kk)=0;
              if jj>3
              if z_pred(kk,5,jj-1)>=z_pred(kk,5,jj-2)  % 这里注意对捕食行为的判断
                  prey_worse(ii,kk)=0;        %是建立在能量不减少，或者能量增加的；
                                              %如果捕食者达到能量下限的情况不变，则会被误认为进食。
              else
                  prey_worse(ii,kk)=1;

              end
              end
           elseif d_shap_3_d<=r_s && d_shap_3_d>=r_c  %视野内有捕食者，但自己没被捕杀
              val = d_shap_3_d;
              chase_12(ii,kk)=val;
              valx = 1/(chase_12(ii,kk)-r_c)^2*...
                     (z_prey(ii,1,jj-1)-shap_3_d(kk,1,jj-1));
              valy = 1/(chase_12(ii,kk)-r_c)^2*...
                     (z_prey(ii,2,jj-1)-shap_3_d(kk,2,jj-1));
              datax(ii,kk)=valx;
              datay(ii,kk)=valy;  
              prey_normal(ii,kk)=0;
              if jj>3
              if z_pred(kk,5,jj-1)>=z_pred(kk,5,jj-2)  % 这里注意对捕食行为的判断
                  prey_worse(ii,kk)=0;        %是建立在能量不减少，或者能量增加的；
                                              %如果捕食者达到能量下限的情况不变，则会被误认为进食。
              else
                  prey_worse(ii,kk)=1;

              end
              end
           elseif d_shap_4_d<=r_s && d_shap_4_d>=r_c  %视野内有捕食者，但自己没被捕杀
              val = d_shap_4_d;
              chase_12(ii,kk)=val;
              valx = 1/(chase_12(ii,kk)-r_c)^2*...
                     (z_prey(ii,1,jj-1)-shap_4_d(kk,1,jj-1));
              valy = 1/(chase_12(ii,kk)-r_c)^2*...
                     (z_prey(ii,2,jj-1)-shap_4_d(kk,2,jj-1));
              datax(ii,kk)=valx;
              datay(ii,kk)=valy;  
              prey_normal(ii,kk)=0;
              if jj>3
              if z_pred(kk,5,jj-1)>=z_pred(kk,5,jj-2)  % 这里注意对捕食行为的判断
                  prey_worse(ii,kk)=0;        %是建立在能量不减少，或者能量增加的；
                                              %如果捕食者达到能量下限的情况不变，则会被误认为进食。
              else
                  prey_worse(ii,kk)=1;

              end
              end
          elseif d_12<r_c              % 距离已经小于被捕杀距离
              z_prey(ii,5,jj)=-1000;  %猎物能量设置为-1000来表示被捕杀
              chase_12(ii,kk)=1000;
              valx = 1/(chase_12(ii,kk)-r_c)^2*...
                     (z_prey(ii,1,jj-1)-z_pred(kk,1,jj-1));
              valy = 1/(chase_12(ii,kk)-r_c)^2*...
                     (z_prey(ii,2,jj-1)-z_pred(kk,2,jj-1));
              datax(ii,kk)=valx;
              datay(ii,kk)=valy;
              prey_normal(ii,kk)=0;
              prey_worse(ii,kk)=0; 
          elseif d_shap_1_d<r_c              % 距离已经小于被捕杀距离
              z_prey(ii,5,jj)=-1000;  %猎物能量设置为-1000来表示被捕杀
              chase_12(ii,kk)=1000;
              valx = 1/(chase_12(ii,kk)-r_c)^2*...
                     (z_prey(ii,1,jj-1)-shap_1_d(kk,1,jj-1));
              valy = 1/(chase_12(ii,kk)-r_c)^2*...
                     (z_prey(ii,2,jj-1)-shap_1_d(kk,2,jj-1));
              datax(ii,kk)=valx;
              datay(ii,kk)=valy;
              prey_normal(ii,kk)=0;
              prey_worse(ii,kk)=0; 
          elseif d_shap_2_d<r_c              % 距离已经小于被捕杀距离
              z_prey(ii,5,jj)=-1000;  %猎物能量设置为-1000来表示被捕杀
              chase_12(ii,kk)=1000;
              valx = 1/(chase_12(ii,kk)-r_c)^2*...
                     (z_prey(ii,1,jj-1)-shap_2_d(kk,1,jj-1));
              valy = 1/(chase_12(ii,kk)-r_c)^2*...
                     (z_prey(ii,2,jj-1)-shap_2_d(kk,2,jj-1));
              datax(ii,kk)=valx;
              datay(ii,kk)=valy;
              prey_normal(ii,kk)=0;
              prey_worse(ii,kk)=0; 
          elseif d_shap_3_d<r_c              % 距离已经小于被捕杀距离
              z_prey(ii,5,jj)=-1000;  %猎物能量设置为-1000来表示被捕杀
              chase_12(ii,kk)=1000;
              valx = 1/(chase_12(ii,kk)-r_c)^2*...
                     (z_prey(ii,1,jj-1)-shap_3_d(kk,1,jj-1));
              valy = 1/(chase_12(ii,kk)-r_c)^2*...
                     (z_prey(ii,2,jj-1)-shap_3_d(kk,2,jj-1));
              datax(ii,kk)=valx;
              datay(ii,kk)=valy;
              prey_normal(ii,kk)=0;
              prey_worse(ii,kk)=0; 
          elseif d_shap_4_d<r_c              % 距离已经小于被捕杀距离
              z_prey(ii,5,jj)=-1000;  %猎物能量设置为-1000来表示被捕杀
              chase_12(ii,kk)=1000;
              valx = 1/(chase_12(ii,kk)-r_c)^2*...
                     (z_prey(ii,1,jj-1)-shap_4_d(kk,1,jj-1));
              valy = 1/(chase_12(ii,kk)-r_c)^2*...
                     (z_prey(ii,2,jj-1)-shap_4_d(kk,2,jj-1));
              datax(ii,kk)=valx;
              datay(ii,kk)=valy;
              prey_normal(ii,kk)=0;
              prey_worse(ii,kk)=0; 
          else % 视野内没有捕食者时，猎物静止不动，静止不动时，不会产生能量消耗
              chase_12(ii,kk)=1000;
              valx = 1/(chase_12(ii,kk)-r_c)^2;
              valy = 1/(chase_12(ii,kk)-r_c)^2;
              datax(ii,kk)=valx;
              datay(ii,kk)=valy;
              prey_normal(ii,kk)=1;
              prey_worse(ii,kk)=0;
          end
    end
              ctx(ii)=sum(datax(ii,:));
              cty(ii)=sum(datay(ii,:));
              f_ct_x(ii)=ctx(ii)/((ctx(ii))^2+(cty(ii))^2);  % 否则，对附近猎物按距离加权求和
              f_ct_y(ii)=cty(ii)/((ctx(ii))^2+(cty(ii))^2);
      %% 判断猎物所处情形
              normal(ii)=all(prey_normal(ii,:)==1); 
              worse(ii)=any(prey_worse(ii,:)==1); 
              situ1 = [normal(ii),worse(ii)];
              better(ii)=all(situ1==0);
              % 猎物视野内没有捕食者，或所有捕食者正在进食，则静止不动
              
              situ2 = [normal(ii),better(ii)];
              prey_stop(ii) = 1-any(situ2==1);        % 捕食者进食时，猎物不动
%               prey_stop(ii) = 1 -normal(ii);        %捕食者进食时，猎物逃跑
     

          uix_1(ii)=beta*pos_x_11(ii)+alfa*spe_x_11(ii)+gamma*f_ct_x(ii);
          uiy_1(ii)=beta*pos_y_11(ii)+alfa*spe_y_11(ii)+gamma*f_ct_y(ii);
      
          v_t(ii) = v_to*sigmf(z_prey(ii,5,jj-1),[10 0.5]);  %速度更新方程
          
          pix_1(ii)=v_t(ii)*uix_1(ii)/(sqrt(uix_1(ii)^2+uiy_1(ii)^2));
          piy_1(ii)=v_t(ii)*uiy_1(ii)/(sqrt(uix_1(ii)^2+uiy_1(ii)^2));
          
          % 给速度增加噪音
          xita_prey = 2*noise_prey*pi*rand()-noise_prey*pi;
          pix_1(ii)=pix_1(ii)*cos(xita_prey)-piy_1(ii)*sin(xita_prey);
          piy_1(ii)=piy_1(ii)*cos(xita_prey)+pix_1(ii)*sin(xita_prey);
         
          gix_1(ii)=gix_1(ii)+pix_1(ii)*delta*prey_stop(ii);  % prey_stop，设置猎物的运动使能
          giy_1(ii)=giy_1(ii)+piy_1(ii)*delta*prey_stop(ii);
%% 猎物的能量更新程序
dot_prey(ii) = z_prey(ii,3,jj-1)*pix_1(ii)+z_prey(ii,4,jj-1)*piy_1(ii);
module_prey(ii) = sqrt((z_prey(ii,3,jj-1))^2+(z_prey(ii,4,jj-1))^2)*...
    sqrt((pix_1(ii))^2+(piy_1(ii))^2);
yuxuan_prey(ii)=dot_prey(ii)/module_prey(ii);
theta_prey(ii) = acos(yuxuan_prey(ii));


if z_prey(ii,5,jj)==-1000
    z_prey(ii,5,jj)=-1000;
else 
    z_prey(ii,5,jj) = z_prey(ii,5,jj-1)+(k_1*normal(ii)+...
        k_2*better(ii)+k_3*worse(ii)*((v_t(ii)/v_to)^2+(theta_prey(ii)/pi)^2))*unit_e;  
end  

%% 猎物能量设置上下限
if z_prey(ii,5,jj)<0 && z_prey(ii,5,jj)>-1000
    z_prey(ii,5,jj)=0;
elseif z_prey(ii,5,jj)>max_e
    z_prey(ii,5,jj)=max_e;
elseif z_prey(ii,5,jj)<min_e && z_prey(ii,5,jj)>0
    z_prey(ii,5,jj)=min_e;
end
%% 猎物在自用空间移动    
%     z_prey(ii,1,jj)=gix_1(ii);
%     z_prey(ii,2,jj)=giy_1(ii);
%     z_prey(ii,3,jj)=pix_1(ii);
%     z_prey(ii,4,jj)=piy_1(ii);
%     
    %% 猎物的活动范围设置为周期环境
    if abs(gix_1(ii))>L/2
    gix_1(ii)=-sign(gix_1(ii))*L+gix_1(ii);
    end
    if abs(giy_1(ii))>L/2
    giy_1(ii)=-sign(giy_1(ii))*L+giy_1(ii);
    end
    
    z_prey(ii,1,jj)=gix_1(ii);
    z_prey(ii,2,jj)=giy_1(ii);
    z_prey(ii,3,jj)=pix_1(ii);
    z_prey(ii,4,jj)=piy_1(ii);
    

 
    
        end
        
 %% 猎物的生存景观图       
    if  z_prey(ii,5,jj)==-1000 && z_prey(ii,5,jj)<z_prey(ii,5,jj-1)
        prey_life(ii)=jj;
    elseif z_prey(ii,5,jj)~=-1000 && jj==timestep
        prey_life(ii)=timestep;
    end
    
    end
   
    %% 指标统计
    % 群体存活状态
    num_dead_prey(jj) = length(find(z_prey(:,5,jj)==-1000));%死亡个数
    num_live_prey(jj) = n_1-num_dead_prey(jj); %幸存个数
    
    
    
    % 捕食者的能量变化图：
    total_energy_pred(jj) = sum(z_pred(:,5,jj));
    % 猎物的能量变化图：
    energy_prey =z_prey(:,5,jj);
    energy_prey(energy_prey==-1000)=0;
    total_energy_prey(jj)= sum(energy_prey);
    % 自然能量总和
    nature_energy(jj)=total_energy_prey(jj)+total_energy_pred(jj);
    
    
    
    
    
    
   fprintf('times: %d \n',jj);
end

%===============绘制初始位置和最终位置图====================================

%    plot(z_pred(:,1,1),z_pred(:,2,1),'og');
%    hold on
%    plot(z_pred(:,1,end),z_pred(:,2,end),'og','markerfacecolor','g');
%    hold on
%    
%    plot(z_prey(:,1,1),z_prey(:,2,1),'^r');
%    hold on
%    plot(z_prey(:,1,end),z_prey(:,2,end),'^r','markerfacecolor','r');
%    hold on
 %========================================================================= 
 
 
%===============绘制初始位置和最终位置图====================================
%起始位置
%    plot(z(1:8,1,1),z(1:8,2,1),'og');
%    hold on 
%    plot(z(10:11,1,1),z(10:11,2,1),'og');
%    hold on 
%    plot(z(13:20,1,1),z(13:20,2,1),'og');
%    hold on 
%    plot(z(9,1,1),z(9,2,1),'*r');
%    hold on 
%    plot(z(12,1,1),z(12,2,1),'*r');
%    hold on 
% %最终位置
%    plot(z(9,1,end),z(9,2,end),'*r','markerfacecolor','r');
%    hold on
%    plot(z(12,1,end),z(12,2,end),'*r','markerfacecolor','r');
%    hold on
%    plot(z(1:8,1,end),z(1:8,2,end),'og','markerfacecolor','g');
%    hold on 
%    plot(z(10:11,1,end),z(10:11,2,end),'og','markerfacecolor','g');
%    hold on 
%    plot(z(13:20,1,end),z(13:20,2,end),'og','markerfacecolor','g');
%    hold on 
   
 %% 计算相应评价指标  
 
 num_dead = length(find(z_prey(:,5,end)==-1000));%死亡个数
 ratio_live =1- num_dead/n_1             %存活率
 
 max_prey_life = max(prey_life);   %种群生存时间
 min_prey_life = min(prey_life) ;  %捕捉到第一个猎物的时间
 ave_prey_life = mean(prey_life)    %个体的平均生存时间
 
 

 total_energy_prey_left=total_energy_prey(timestep);  %猎物种群能量剩余
 ave_energy_prey_left = total_energy_prey_left/(n_1-num_dead)
 total_energy_pred_left=total_energy_pred(timestep);  %捕食者种群能量剩余
 ave_energy_pred_left = total_energy_pred_left/n_2
 
 nature_energy_left = total_energy_prey_left+total_energy_pred_left %自然能量剩余
 nature_energy_waste = nature_energy(1) - nature_energy(timestep); %捕食对抗的能量浪费
 %% 能量变化图
 figure(1)
 plot(total_energy_prey);
 title('猎物总能量变化曲线')
 
 figure(2)
 plot(total_energy_pred);
 title('捕食者总能量变化曲线')
 
 figure(3)
 plot(nature_energy);
 title('自然总能量变化曲线')
 
 figure(4)
 plot(num_live_prey);
 title('猎物生存个数')
 xlabel('time step');
 ylabel('alive prey');
%  axis equal
%  xlabel('x/m');
%  ylabel('y/m');


% 绘制轨迹曲线=======================================================

% aaa=z_pred(:,1,:);
%  bbb=z_pred(:,2,:);
%  aaa=reshape(aaa,n_2,timestep);
%  bbb=reshape(bbb,n_2,timestep);
%  for kkk=1:n_2
%   plot(aaa(kkk,:),bbb(kkk,:),':b');
%   hold on
%  end
%  
%  ccc=z_prey(:,1,:);
%  ddd=z_prey(:,2,:);
%  ccc=reshape(ccc,n_1,timestep);
%  ddd=reshape(ddd,n_1,timestep);
%  for kkk=1:n_1
%   plot(ccc(kkk,:),ddd(kkk,:),':k');
%   hold on
%  end

 