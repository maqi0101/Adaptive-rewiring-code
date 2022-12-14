function dydt = LV_PHM(t,y,Sp,Sh,Sm,r,beta_P,beta_M,beta_H,alpha_PH,gamma_PM,e,h)
dydt=zeros(Sp+Sh+Sm,1);%Sp:plant;Sh:herbivore;Sm:bee，列向量
Pi=y(1:Sp);% plant
Hi=y(Sp+1:Sp+Sh);%herbivore
Mi=y(Sp+Sh+1:Sp+Sh+Sm);%pollinator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
matrix_PH=logical(alpha_PH);
matrix_PM=logical(gamma_PM);
%食草动物和植物间相互作用
F_ph=matrix_PH.*(Pi*Hi')./(1+h*repmat((matrix_PH'*Pi)',[Sp,1]));
%%
%传粉者从植物获得效应
E1_pm=matrix_PM.*(Pi*Mi')./(1+h*repmat((matrix_PM'*Pi)',[Sp,1]));
%%
%植物从传粉者获得效应
E2_mp=matrix_PM'.*(Mi*Pi')./(1+h*repmat((matrix_PM*Mi)',[Sm,1]));
%%
r_P=r(1:Sp);
r_H=r(Sp+1:Sp+Sh);
r_M=r(Sp+Sh+1:Sp+Sh+Sm);
%%
dPi=Pi.*(r_P-beta_P*Pi)+(sum(gamma_PM'.*E2_mp))'-sum(alpha_PH.*F_ph,2);
dHi=Hi.*(r_H-beta_H*Hi)+(sum(e*alpha_PH.*F_ph))';
dMi=Mi.*(r_M-beta_M*Mi)+(sum(gamma_PM.*E1_pm))';
dydt(1:Sp)=dPi;
dydt(Sp+1:Sp+Sh)=dHi;
dydt(Sp+Sh+1:Sp+Sh+Sm)=dMi;