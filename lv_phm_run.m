clc;
clear; 
Sp=30;%%%%%%plant的物种数
Sh=30;%%%%%%%食草动物物种数
Sm=30;%%%%%%%%蜜蜂物种数
mu_r=1;%%用来调节物种自然增长率
r=mu_r*ones(Sp+Sh+Sm,1); %%物种自然增长率
h=0.1;%%%%功能反应的半饱和常数
e=0.8;%%%%%植物到捕食者的转化系数
%%
eta=1;   %%用来调节断开链接概率的参数
trSpan = 1;  %%%%生态位0-1
nw=0.1;    %%%生态位宽度
int_comp=0.1;   %%%竞争强度
int_mut=0.1;   %%%互惠强度
int_ant=0.1;   %%%捕食强度
%%
trait=rand(Sp+Sh+Sm,1);%%%每个物种的生态位的中心位置
Cmut = 0.15;   %互惠连接度
Cant = 0.15;   %捕食连接度
%%
beta_P = zeros(Sp);  %%植物之间的竞争
beta_H = zeros(Sh); %%食草动物之间的竞争
beta_M = zeros(Sm); %%传粉者之间的竞争
%
for i=1:Sp
    for j=i+1:Sp        
        beta_P(i,j)=int_comp*overlap(trait(i),trait(j),nw,trSpan);%%植物的种间竞争
        beta_P(j,i)=beta_P(i,j);        
    end
end
beta_P = beta_P + eye(Sp);%%植物的种内竞争
%
for i=1:Sh
    for j=i+1:Sh        
        beta_H(i,j)=int_comp*overlap(trait(Sp+i),trait(Sp+j),nw,trSpan); %%食草动物的种间竞争  
        beta_H(j,i)=beta_H(i,j);        
    end
end
beta_H=beta_H + eye(Sh);%%食草动物的种内竞争
%
for i=1:Sm
    for j=i+1:Sm        
        beta_M(i,j)=int_comp*overlap(trait(Sp+Sh+i),trait(Sp+Sh+j),nw,trSpan); %%传粉者的种间竞争  
        beta_M(j,i)=beta_M(i,j);        
    end
end
beta_M=beta_M + eye(Sm);%%传粉者的种内竞争
%%
alpha_PH=zeros(Sp,Sh);
for i=1:Sp
    for j=1:Sh
        if rand < Cant
            alpha_PH(i,j)=int_ant*overlap(trait(i),trait(Sp+j),nw,trSpan);%%捕食的相互作用
        end
    end
end
%
gamma_PM=zeros(Sp,Sm);
for i=1:Sp
    for j=1:Sm  
        if rand < Cmut
            gamma_PM(i,j)=int_mut*overlap(trait(i),trait(Sp+Sh+j),nw,trSpan);  %%互惠的相互作用
        end
    end
end
%%
counter_max = 100000;%%%%适应性重连的总次数
counter=0;%%%用来记录重连的次数
idstart = 0;%%用来储存每次重连物种的位置
Type = 0;%%%用来选择重连物种所在的行会,1:H->P; 2:M->P
mu_N0 = 0.1;%%用来调节物种初始生物量
N0 = mu_N0*ones(Sp+Sh+Sm,1);  %%%物种初始生物量
to = 0;   %初始时间t0
tf = 20;  % integration time within each time step %%时间间隔
%%
while counter<counter_max
    counter = counter+1;   
    % 挑选一个物种让其进行适应性重连
    Type=randi(2);   %%% 选择重连物种所在的行会,1:H->P; 2:M->P             
    [alpha_PH,gamma_PM,idstart,idend_old,idend_new,gain_old_rw,int_old]=rewiring(Type,N0,Sp,Sh,Sm,alpha_PH,gamma_PM,trait,int_ant,int_mut,nw,trSpan,eta);
    % integration
    [t y] = ode45(@(t,y)LV_PHM(t,y,Sp,Sh,Sm,r,beta_P,beta_M,beta_H,alpha_PH,gamma_PM,e,h),[to to+tf],N0);%%%重连后，运行到稳定
    to = t(end);
    N0=y(end,:)';
    %
    jacob_mat=get_jacmat(N0,Sp,Sh,Sm,r,beta_P,beta_M,beta_H,alpha_PH,gamma_PM,e,h); %雅克比矩阵
    lambda=-max(real(eig(jacob_mat))); % 稳定性（雅克比矩阵特征值负的最大实部）
    if counter>1
        if Type==1
            gain_new_rw=N0(idstart+Sp); %%该物种重连后的生物量
        else
            gain_new_rw=N0(idstart+Sp+Sh); 
        end
        if gain_old_rw > gain_new_rw     % 比较重连前后的生物量，判断是否恢复原来的链接
            [alpha_PH,gamma_PM]=wire_back(alpha_PH,gamma_PM,int_old,idstart,idend_new,idend_old,Type);      
        end
    end
    %
    [nodf,qb,Nm] = cal_structure(alpha_PH); % 计算捕食子网络的嵌套性和模块化
    N_ant=nodf;  %捕食子网络的嵌套性
    Q_ant=qb;    %捕食子网络的模块化
    %
    [nodf,qb,Nm] = cal_structure(gamma_PM); % 计算互惠子网络的嵌套性和模块化
    N_mut=nodf;   %互惠子网络的嵌套性
    Q_mut=qb;     %互惠子网络的模块化
end
