clc;
clear; 
Sp=30;%%%%%%plant��������
Sh=30;%%%%%%%ʳ�ݶ���������
Sm=30;%%%%%%%%�۷�������
mu_r=1;%%��������������Ȼ������
r=mu_r*ones(Sp+Sh+Sm,1); %%������Ȼ������
h=0.1;%%%%���ܷ�Ӧ�İ뱥�ͳ���
e=0.8;%%%%%ֲ�ﵽ��ʳ�ߵ�ת��ϵ��
%%
eta=1;   %%�������ڶϿ����Ӹ��ʵĲ���
trSpan = 1;  %%%%��̬λ0-1
nw=0.1;    %%%��̬λ���
int_comp=0.1;   %%%����ǿ��
int_mut=0.1;   %%%����ǿ��
int_ant=0.1;   %%%��ʳǿ��
%%
trait=rand(Sp+Sh+Sm,1);%%%ÿ�����ֵ���̬λ������λ��
Cmut = 0.15;   %�������Ӷ�
Cant = 0.15;   %��ʳ���Ӷ�
%%
beta_P = zeros(Sp);  %%ֲ��֮��ľ���
beta_H = zeros(Sh); %%ʳ�ݶ���֮��ľ���
beta_M = zeros(Sm); %%������֮��ľ���
%
for i=1:Sp
    for j=i+1:Sp        
        beta_P(i,j)=int_comp*overlap(trait(i),trait(j),nw,trSpan);%%ֲ����ּ侺��
        beta_P(j,i)=beta_P(i,j);        
    end
end
beta_P = beta_P + eye(Sp);%%ֲ������ھ���
%
for i=1:Sh
    for j=i+1:Sh        
        beta_H(i,j)=int_comp*overlap(trait(Sp+i),trait(Sp+j),nw,trSpan); %%ʳ�ݶ�����ּ侺��  
        beta_H(j,i)=beta_H(i,j);        
    end
end
beta_H=beta_H + eye(Sh);%%ʳ�ݶ�������ھ���
%
for i=1:Sm
    for j=i+1:Sm        
        beta_M(i,j)=int_comp*overlap(trait(Sp+Sh+i),trait(Sp+Sh+j),nw,trSpan); %%�����ߵ��ּ侺��  
        beta_M(j,i)=beta_M(i,j);        
    end
end
beta_M=beta_M + eye(Sm);%%�����ߵ����ھ���
%%
alpha_PH=zeros(Sp,Sh);
for i=1:Sp
    for j=1:Sh
        if rand < Cant
            alpha_PH(i,j)=int_ant*overlap(trait(i),trait(Sp+j),nw,trSpan);%%��ʳ���໥����
        end
    end
end
%
gamma_PM=zeros(Sp,Sm);
for i=1:Sp
    for j=1:Sm  
        if rand < Cmut
            gamma_PM(i,j)=int_mut*overlap(trait(i),trait(Sp+Sh+j),nw,trSpan);  %%���ݵ��໥����
        end
    end
end
%%
counter_max = 100000;%%%%��Ӧ���������ܴ���
counter=0;%%%������¼�����Ĵ���
idstart = 0;%%��������ÿ���������ֵ�λ��
Type = 0;%%%����ѡ�������������ڵ��л�,1:H->P; 2:M->P
mu_N0 = 0.1;%%�����������ֳ�ʼ������
N0 = mu_N0*ones(Sp+Sh+Sm,1);  %%%���ֳ�ʼ������
to = 0;   %��ʼʱ��t0
tf = 20;  % integration time within each time step %%ʱ����
%%
while counter<counter_max
    counter = counter+1;   
    % ��ѡһ���������������Ӧ������
    Type=randi(2);   %%% ѡ�������������ڵ��л�,1:H->P; 2:M->P             
    [alpha_PH,gamma_PM,idstart,idend_old,idend_new,gain_old_rw,int_old]=rewiring(Type,N0,Sp,Sh,Sm,alpha_PH,gamma_PM,trait,int_ant,int_mut,nw,trSpan,eta);
    % integration
    [t y] = ode45(@(t,y)LV_PHM(t,y,Sp,Sh,Sm,r,beta_P,beta_M,beta_H,alpha_PH,gamma_PM,e,h),[to to+tf],N0);%%%���������е��ȶ�
    to = t(end);
    N0=y(end,:)';
    %
    jacob_mat=get_jacmat(N0,Sp,Sh,Sm,r,beta_P,beta_M,beta_H,alpha_PH,gamma_PM,e,h); %�ſ˱Ⱦ���
    lambda=-max(real(eig(jacob_mat))); % �ȶ��ԣ��ſ˱Ⱦ�������ֵ�������ʵ����
    if counter>1
        if Type==1
            gain_new_rw=N0(idstart+Sp); %%�������������������
        else
            gain_new_rw=N0(idstart+Sp+Sh); 
        end
        if gain_old_rw > gain_new_rw     % �Ƚ�����ǰ������������ж��Ƿ�ָ�ԭ��������
            [alpha_PH,gamma_PM]=wire_back(alpha_PH,gamma_PM,int_old,idstart,idend_new,idend_old,Type);      
        end
    end
    %
    [nodf,qb,Nm] = cal_structure(alpha_PH); % ���㲶ʳ�������Ƕ���Ժ�ģ�黯
    N_ant=nodf;  %��ʳ�������Ƕ����
    Q_ant=qb;    %��ʳ�������ģ�黯
    %
    [nodf,qb,Nm] = cal_structure(gamma_PM); % ���㻥���������Ƕ���Ժ�ģ�黯
    N_mut=nodf;   %�����������Ƕ����
    Q_mut=qb;     %�����������ģ�黯
end
