function [alpha_PH,gamma_PM,idstart,idend_old,idend_new,gain_old,int_old]=rewiring(Type,N,Sp,Sh,Sm,alpha_PH,gamma_PM,trait,int_ant,int_mut,nw,trSpan,eta)

% rewiring one link

if (Type==1)  % rewiring H->P %%%从食草动物行会中确定
    idH=randi(Sh);%%%随机从食草动物确定一个物种
    while isempty(find(alpha_PH(:,idH),1))  %%判断该动物是否有连接的植物伙伴，没有的话重新选择一个动物
        idH=randi(Sh);
    end
    
    gain_old=N(Sp+idH);  %记录该动物重连前的生物量
    
    idP_nz=find(alpha_PH(:,idH));    %找出该动物的全部植物伙伴
    idP_old=idP_nz(randi(length(idP_nz)));%%从中随机选择一个植物
    
    int_old=alpha_PH(idP_old,idH);%%%记录他们的作用强度
    
    idP_new=idP_old;
    deg_old=nnz(alpha_PH(idP_old,:));%%%找出该植物一共有多少连接伙伴
    pr=1/deg_old^eta;  % rewiring probability %%%断开的概率
    if rand>pr   
        idP_z=find(alpha_PH(:,idH)==0);%%找出该动物的未连接的植物
        if ~isempty(idP_z) 
        idP_new=idP_z(randi(length(idP_z))); %%从未连接的植物中随机选择一个植物
        % rewire to a new species
        alpha_PH(idP_old,idH)=0;%%断开原来的连接
        
        alpha_PH(idP_new,idH)=int_ant*overlap(trait(idP_new),trait(Sp+idH),nw,trSpan); %%%把新的植物与该动物连接，并更新他们的作用强度
        end
    end    
    
    idstart=idH;%%%记录该动物的编号
    idend_old=idP_old;%%记录旧连接植物的编号
    idend_new=idP_new;%%%记录新连接植物的编号
       
else     % rewiring M->P %%%从传粉者行会中确定
    idM=randi(Sm);%%%随机从传粉者确定一个物种
    while isempty(find(gamma_PM(:,idM),1))  % %%判断该动物是否有连接的植物伙伴，没有的话重新选择一个动物
        idM=randi(Sm);
    end
    
    gain_old=N(Sp+Sh+idM); %记录该动物重连前的生物量
    
    idP_nz=find(gamma_PM(:,idM));%找出该动物的全部植物伙伴
    idP_old=idP_nz(randi(length(idP_nz)));%%从中随机选择一个植物
    
    int_old=gamma_PM(idP_old,idM);%%%记录他们的作用强度
    
    idP_new=idP_old;
    deg_old=nnz(gamma_PM(idP_old,:));%%%找出该植物一共有多少连接伙伴
    pr=1/deg_old^eta;%%%断开的概率
    if rand>pr   
        idP_z=find(gamma_PM(:,idM)==0);%%找出该动物的未连接的植物
        if ~isempty(idP_z)        
        idP_new=idP_z(randi(length(idP_z)));%%从未连接的植物中随机选择一个植物
        % rewire to a new species
        gamma_PM(idP_old,idM)=0; %%断开原来的连接
        gamma_PM(idP_new,idM)=int_mut*overlap(trait(idP_new),trait(Sp+Sh+idM),nw,trSpan);%%%把新的植物与该动物连接，并更新他们的作用强度
        end
    end
    
    idstart=idM;
    idend_old=idP_old;
    idend_new=idP_new;
    
end

