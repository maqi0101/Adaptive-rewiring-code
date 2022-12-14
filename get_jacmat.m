function jac_mat = get_jacmat(ymea,Sp,Sh,Sm,r,beta_P,beta_M,beta_H,alpha_PH,gamma_PM,e,h)
%%
PP=ymea(1:Sp);% plant
HH=ymea(Sp+1:Sp+Sh);%herbivore
MM=ymea(Sp+Sh+1:Sp+Sh+Sm);%pollinator
r_P=r(1:Sp);
r_H=r(Sp+1:Sp+Sh);
r_M=r(Sp+Sh+1:Sp+Sh+Sm);
mat_PH=logical(alpha_PH);
mat_PM=logical(gamma_PM);
mat_MP=mat_PM';
mat_HP=mat_PH';
%%
dif_PP=zeros(Sp,Sp);
for i=1:Sp
    for j=1:Sp
        dif_PP(i,j)=PP(i)*(-beta_P(i,j)+h*alpha_PH(i,:).*mat_PH(j,:)./((1+h*mat_HP*PP)'.^2)*HH);
    end
end
%%
dif_PH=zeros(Sp,Sh);
for i=1:Sp
    for j=1:Sh
        dif_PH(i,j)=-PP(i)*alpha_PH(i,j)/(1+h*mat_HP(j,:)*PP);
    end
end
%%
dif_PM=zeros(Sp,Sm);
for i=1:Sp
    for j=1:Sm
        dif_PM(i,j)=PP(i)/(1+h*mat_PM(i,:)*MM)*(gamma_PM(i,j)+mat_PM(i,j)*h*(r_P(i)-beta_P(i,:)*PP-alpha_PH(i,:)./((1+h*mat_HP*PP)')*HH));
    end
end
%%
dif_HP=zeros(Sh,Sp);
for i=1:Sh
    for j=1:Sp
        dif_HP(i,j)=HH(i)*(e*alpha_PH(j,i)+mat_HP(i,j)*h*(r_H(i)-beta_H(i,:)*HH))/(1+h*mat_HP(i,:)*PP);
    end
end
%%
dif_HH=zeros(Sh,Sh);
for i=1:Sh
    for j=1:Sh
        dif_HH(i,j)=-HH(i)*beta_H(i,j);
    end
end
%%
dif_HM=zeros(Sh,Sm);
%%
dif_MP=zeros(Sm,Sp);
for i=1:Sm
    for j=1:Sp
        dif_MP(i,j)=MM(i)*(gamma_PM(j,i)+mat_MP(i,j)*h*(r_M(i)-beta_M(i,:)*MM))/(1+h*mat_MP(i,:)*PP);
    end
end
%%
dif_MH=zeros(Sm,Sh);
%%
dif_MM=zeros(Sm,Sm);
for i=1:Sm
    for j=1:Sm
        dif_MM(i,j)=-MM(i)*beta_M(i,j);
    end
end
%%
jac_mat=[dif_PP,dif_PH,dif_PM;
         dif_HP,dif_HH,dif_HM;
         dif_MP,dif_MH,dif_MM];