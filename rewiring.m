function [alpha_PH,gamma_PM,idstart,idend_old,idend_new,gain_old,int_old]=rewiring(Type,N,Sp,Sh,Sm,alpha_PH,gamma_PM,trait,int_ant,int_mut,nw,trSpan,eta)

% rewiring one link

if (Type==1)  % rewiring H->P %%%��ʳ�ݶ����л���ȷ��
    idH=randi(Sh);%%%�����ʳ�ݶ���ȷ��һ������
    while isempty(find(alpha_PH(:,idH),1))  %%�жϸö����Ƿ������ӵ�ֲ���飬û�еĻ�����ѡ��һ������
        idH=randi(Sh);
    end
    
    gain_old=N(Sp+idH);  %��¼�ö�������ǰ��������
    
    idP_nz=find(alpha_PH(:,idH));    %�ҳ��ö����ȫ��ֲ����
    idP_old=idP_nz(randi(length(idP_nz)));%%�������ѡ��һ��ֲ��
    
    int_old=alpha_PH(idP_old,idH);%%%��¼���ǵ�����ǿ��
    
    idP_new=idP_old;
    deg_old=nnz(alpha_PH(idP_old,:));%%%�ҳ���ֲ��һ���ж������ӻ��
    pr=1/deg_old^eta;  % rewiring probability %%%�Ͽ��ĸ���
    if rand>pr   
        idP_z=find(alpha_PH(:,idH)==0);%%�ҳ��ö����δ���ӵ�ֲ��
        if ~isempty(idP_z) 
        idP_new=idP_z(randi(length(idP_z))); %%��δ���ӵ�ֲ�������ѡ��һ��ֲ��
        % rewire to a new species
        alpha_PH(idP_old,idH)=0;%%�Ͽ�ԭ��������
        
        alpha_PH(idP_new,idH)=int_ant*overlap(trait(idP_new),trait(Sp+idH),nw,trSpan); %%%���µ�ֲ����ö������ӣ����������ǵ�����ǿ��
        end
    end    
    
    idstart=idH;%%%��¼�ö���ı��
    idend_old=idP_old;%%��¼������ֲ��ı��
    idend_new=idP_new;%%%��¼������ֲ��ı��
       
else     % rewiring M->P %%%�Ӵ������л���ȷ��
    idM=randi(Sm);%%%����Ӵ�����ȷ��һ������
    while isempty(find(gamma_PM(:,idM),1))  % %%�жϸö����Ƿ������ӵ�ֲ���飬û�еĻ�����ѡ��һ������
        idM=randi(Sm);
    end
    
    gain_old=N(Sp+Sh+idM); %��¼�ö�������ǰ��������
    
    idP_nz=find(gamma_PM(:,idM));%�ҳ��ö����ȫ��ֲ����
    idP_old=idP_nz(randi(length(idP_nz)));%%�������ѡ��һ��ֲ��
    
    int_old=gamma_PM(idP_old,idM);%%%��¼���ǵ�����ǿ��
    
    idP_new=idP_old;
    deg_old=nnz(gamma_PM(idP_old,:));%%%�ҳ���ֲ��һ���ж������ӻ��
    pr=1/deg_old^eta;%%%�Ͽ��ĸ���
    if rand>pr   
        idP_z=find(gamma_PM(:,idM)==0);%%�ҳ��ö����δ���ӵ�ֲ��
        if ~isempty(idP_z)        
        idP_new=idP_z(randi(length(idP_z)));%%��δ���ӵ�ֲ�������ѡ��һ��ֲ��
        % rewire to a new species
        gamma_PM(idP_old,idM)=0; %%�Ͽ�ԭ��������
        gamma_PM(idP_new,idM)=int_mut*overlap(trait(idP_new),trait(Sp+Sh+idM),nw,trSpan);%%%���µ�ֲ����ö������ӣ����������ǵ�����ǿ��
        end
    end
    
    idstart=idM;
    idend_old=idP_old;
    idend_new=idP_new;
    
end

