clc; clear all;
testmax=20;
A_itr=zeros(testmax,1);A_time_used=A_itr;B_itr=A_itr;B_time_used=A_itr;
for t=1:testmax
[G,ConstrA]=Case_test(8,1000,[300,300,300]);
[X,z,info,BB_lbfgs_itr,BB_lbfgs_time_used] =  CaliMatLbfgsbb_adapt1(G,ConstrA,[]);
[X1,z1,F, BB_lbfgsssn_itr,BB_lbfgsssn_time_used]= CaliMatLbfgsbb_adapt2(G,ConstrA,[]);

%获取实验结果并准备输出
A_itr(t)        =BB_lbfgs_itr;
A_time_used(t)  =BB_lbfgs_time_used;
B_itr(t)        =BB_lbfgsssn_itr;
B_time_used(t)  =BB_lbfgsssn_time_used;
if t==testmax
    Test_Information=table((1:testmax)', A_itr,A_time_used,B_itr,B_time_used,...
        'VariableNames',{'Test_Number','LBFGS_Iterations',...
        'LBFGS_Time_Used','LBFGSSSN_Iterations','LBFGSSSN_Iime_Used'})
    
     path=['Log\Test_Information_',num2str(t),'.xls'];
     writetable(Test_Information,path);
%      xlswrite(path,Test_Information)
end


if  BB_lbfgsssn_time_used > BB_lbfgs_time_used || BB_lbfgsssn_itr>BB_lbfgs_itr
%         pathOPTIONS=['Log\OPTIONS',t,'.mat'];
        pathConstrA=['Log\ConstrA',num2str(t),'.mat'];
        pathG=['Log\G',num2str(t),'.mat'];
%         save(pathOPTIONS,OPTIONS) ;
        save(pathConstrA,'ConstrA') ;
        save(pathG, 'G') ;
        
        %%输出
        Test_Information=table((1:testmax)', A_itr,A_time_used,B_itr,B_time_used,...
        'VariableNames',{'Test_Number','LBFGS_Iterations',...
        'LBFGS_Time_Used','LBFGSSSN_Iterations','LBFGSSSN_Iime_Used'})
        path=['Log\Test_Information',num2str(t),'.xls'];
        writetable(Test_Information,path);
%       xlswrite(path,Test_Information)
         continue
end
end

%%
% % load ConstrA;
% %
shar_J=(ConstrB.Jl==ConstrB.Ju);
shar_l=ConstrB.l(shar_J);
shar_u=ConstrB.u(shar_J);
plot([shar_l,shar_u],'o');
