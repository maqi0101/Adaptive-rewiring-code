function [alpha_PH,gamma_PM] = wire_back(alpha_PH,gamma_PM,int_old,idstart,idend_new,idend_old,Type)

if Type==1 
    alpha_PH(idend_new,idstart)=0;%%将新的连接给断开  
    alpha_PH(idend_old,idstart)=int_old;  %%恢复原来的连接
else
    gamma_PM(idend_new,idstart)=0; 
    gamma_PM(idend_old,idstart)=int_old; 
end 