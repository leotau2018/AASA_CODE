%% 本文档用来集中运行代码，测试。%% test_BB_LBFGS_SN_numericresult
% 清除缓存。
clear all;
experiments_date = string(year(date))+string(month(date))+string(day(date));
result_dir = ['TestResult',char(experiments_date)];
mkdir(result_dir);
%% 主程序开始，设置测试参数。
%%%%%%%这里设定测试问题的参数
%%%%%%%取消table1 2 3 4相应块的注释，运行即可。
%% tabl 1 %约12个小时的计算时间
% matrix_n=[400,500,600,700,800,900];constraint_n =[100,200,300,400,500]; CASETYPE=[1];
%% tabl 2 %约12个小时的计算时间
% matrix_n=[400,500,600,700,800,900];constraint_n =[100,200,300,400,500]; CASETYPE=[2];
%% tabl 3 %约8个小时的计算时间
% matrix_n=[600,700,800,900];        constraint_n =[500];                 CASETYPE=[5,6,7,8];
%% tabl 4 %约60个小时的计算时间
% matrix_n=[1000,1500,2000];         constraint_n =[300,600,800];         CASETYPE=[5,6,7,8];  %[1,2,]
%% %%%%%%%%%%%%%%%% generate Figure Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% synthetic easy example
% matrix_n=[1000]; constraint_n =[500];   CASETYPE=[5,6];    diary_on = 1;
%% synthetic hard example
% matrix_n=[2000]; constraint_n =[500];        CASETYPE=[5,6];    diary_on = 1;
%% Real Cases
% matrix_n= [792,  1187];     constraint_n = [792, 400];      CASETYPE=[3, 2];
% real_typeG = [1,2];     diary_on = 1;
%%%% [G,ConstrA]=case_real_test(1,3,792);            real_case = 1;       % Real Case E7
%%% [G,ConstrA]=case_real_test(2,2,[400,400,400]);  real_case = 2;       % Real Case E8
%% other cases
% matrix_n=[10,500,1000,1500,2000,2500,3000];   %矩阵的维数
% constraint_n =[3,100,200,300,400,500,600,800,100,1200,1400,1600,1800,2000,2200];
% %每一行约束的个数
% CASETYPE=[1,2];  %测试的问题类别

% matrix_n=[400,500,600,700,800,900];constraint_n =[100,200,250,300,400,500]; CASETYPE=[1,2,5,6,7,8];
% matrix_n=[1000,1200,1400,1600,1800,2000];constraint_n =[300,400,500,600,700,800]; CASETYPE=[5,6,7,8];

% matrix_n=[2000];constraint_n =[800]; CASETYPE=[7,8];
% matrix_n=[10,20];constraint_n =[3,5,10]; CASETYPE=[7,8];
matrix_n=[400,900];constraint_n =[100,500]; CASETYPE=[2];
%-------------------------------------------------------------------------
if exist('diary_on','var') == 1
    DIARY = ['Diary',char(experiments_date)];
    DIARY_log_name = [DIARY,'/log_',num2str(matrix_n),'_',num2str(constraint_n),'_',datestr(now,30)];
    mkdir(DIARY); 
    diary(DIARY_log_name);
    diary on;%建立保存输出到命令窗口的文本的文档,
end
% 如果要记录输出的每一步，可以取消这行注释，可能会稍微降低运行效率
%% Run and record
list_n=0;
for mi =1:length(matrix_n)
    for ci=1:length(constraint_n)
        for typ=1:length(CASETYPE)
            if constraint_n(ci)>=matrix_n(mi) && exist('real_typeG','var')==0 
                continue
            end
                casetype=CASETYPE(typ) ;               %测试问题的类型，1为带状对角矩阵约束；2为约束位置随机分布。
                n=matrix_n(mi);                        %测试问题的维数,测试取值1000，1500，2000
                constr_n=[1,1,1]*constraint_n(ci);     %约束个数，三维非负向量，且小于n。
                
                test_n =0;test_N=1;  %%%%变量，test_n =1:test_N;%参数，test_N表示每个instance上总的测试次数
                
                %test_N==1,每个instance上测试1次，所有案例测试完之后保存在一个excel中
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % tau = .0e-8;
                % OPTIONS.tau = tau;
                OPTIONS.tol = 1e-6;  
                OPTIONS.M   = 3; % OPTIONS.M为L-bfgs中的参数L
                OPTIONS.Cumpt_Condn = 1;
                if n >= 1000
                    OPTIONS.maxtime=3600;   % 大规模时，设置时间上限为3600秒
                else
                    OPTIONS.maxitr = 1500; % 中规模时，设置时间上限为1500秒
                end
                if exist('diary_on','var') == 1, OPTIONS.maxtime=100000; OPTIONS.maxit=100000; end
                if exist('real_typeG','var') 
                    if mi == ci && mi == typ
                        [G,ConstrA]=case_real_test(real_typeG(mi),CASETYPE(typ),constraint_n(ci)); 
                        %continue
                    else
                        continue
                    end
                else
                    [G,ConstrA]=Case_test(casetype, n, constr_n);
                end
                
                
                e = ConstrA.e; I_e = ConstrA.Ie; J_e = ConstrA.Je;
                l = ConstrA.l; I_l = ConstrA.Il; J_l = ConstrA.Jl;
                u = ConstrA.u; I_u = ConstrA.Iu; J_u = ConstrA.Ju;
                k_e=length(I_e);k_l=length(I_l);k_u=length(I_u);
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %比较和测试
                %if n>1000&&constr_n(1)>300&&test_n>2
                if n>=1000 && exist('real_typeG','var')==0 %% 大规模(n>1000)时ADM 和 projbfgs 都不再计算
                    OPTIONS_limit=OPTIONS;OPTIONS_limit.maxtime=2;OPTIONS_limit.maxit=2;%%大规模的时候不再多次计算，太费时间, 特别是 projbfgs（P-BFGS）
                else,OPTIONS_limit=OPTIONS;
                end
                [CADM_X,CADM_Z,CADM_info,CADM_itr,CADM_time_used]=CADM(G,ConstrA,OPTIONS_limit);
                
                [pbfgs_Z,pbfgs_itr, pbfgs_histout,costdata,pbfgs_time_used, pbfgs_flag] = projbfgs(G,ConstrA,OPTIONS_limit);
                
                [BB_lbfgsssn_X,BB_lbfgsssn_z,BB_lbfgsssn_info, BB_lbfgsssn_itr,BB_lbfgsssn_time_used,condnum]...
                    = AASA_BB_Lbfgs_SNN_adapt_Ref(G,ConstrA,OPTIONS);
 
                [newton_X,newton_z,newton_info,newton_itr,f_eval,newton_time_used] = CaliMat1Mex(G,ConstrA,OPTIONS);
                
                [BB_lbfgs_X,BB_lbfgs_z,info,BB_lbfgs_itr,BB_lbfgs_time_used] =  AASA_BB_Lbfgs_SNN_adapt1(G,ConstrA,OPTIONS);
                
                %[Lbfgs2D_X,Lbfgs2D_z,Lbfgs2D_F,Lbfgs2D_itr,Lbfgs2D_time_used]      = CaliMatLbfgs2D(G,ConstrA,OPTIONS);
                
                list_n=list_n+1;
                
                % nrmz=norm([newton_z.e;newton_z.l;newton_z.u]-[BB_lbfgsssn_z.e;BB_lbfgsssn_z.l;BB_lbfgsssn_z.u]);
                [theta0,g_ze,g_zl,g_zu,~] = thetafun(G,e,newton_z.e,I_e,J_e,l,newton_z.l,I_l,J_l,u,newton_z.u,I_u,J_u,n);
                z_lu=[newton_z.l;newton_z.u];
                pg=[-g_ze;proje(z_lu-[g_zl;g_zu],zeros(k_l+k_u,1))-z_lu];
                nrmG0 = norm(pg);
                % nrmG0 = norm([g_ze;g_zl;g_zu], 'fro');
                
                [theta,g_ze,g_zl,g_zu,~] = thetafun(G,e,BB_lbfgs_z.e,I_e,J_e,l,BB_lbfgs_z.l,I_l,J_l,u,BB_lbfgs_z.u,I_u,J_u,n);
                z_lu=[BB_lbfgs_z.l;BB_lbfgs_z.u];
                pg=[-g_ze;proje(z_lu-[g_zl;g_zu],zeros(k_l+k_u,1))-z_lu];
                nrmG = norm(pg);
                % nrmG = norm([g_ze;g_zl;g_zu], 'fro');
                
                [theta1,g_ze,g_zl,g_zu,~] = thetafun(G,e,BB_lbfgsssn_z.e,I_e,J_e,l,BB_lbfgsssn_z.l,I_l,J_l,u,BB_lbfgsssn_z.u,I_u,J_u,n);
                z_lu=[BB_lbfgsssn_z.l;BB_lbfgsssn_z.u];
                pg=[-g_ze;proje(z_lu-[g_zl;g_zu],zeros(k_l+k_u,1))-z_lu];
                nrmG1 = norm(pg);
                % nrmG1 = norm([g_ze;g_zl;g_zu], 'fro');
                
                [theta2,g_ze,g_zl,g_zu,eig_time] = thetafun(G,e,pbfgs_Z(1:k_e),I_e,J_e,l,pbfgs_Z(k_e+1:k_e+k_l),I_l,J_l,u,pbfgs_Z(k_e+k_l+1:k_e+k_l+k_u),I_u,J_u,n);
                z_lu=pbfgs_Z(k_e+1:end);
                pg=[-g_ze;proje(z_lu-[g_zl;g_zu],zeros(k_l+k_u,1))-z_lu];
                nrmG2 = norm(pg);
                
                diftheta=theta0-theta1;
                difnrm=nrmG0-nrmG1;
                
                Newton_itr(list_n)            =newton_itr;
                Newton_time_used(list_n)      =newton_time_used;
                Newton_normG(list_n)          =nrmG0;
                A_Lbfgs_itr(list_n)           =BB_lbfgs_itr;
                A_Lbfgs_time_used(list_n)     =BB_lbfgs_time_used;
                A_Lbfgs_normG(list_n)         =nrmG;
                A_ada_itr(list_n)             =BB_lbfgsssn_itr.total_itr;
                A_ada_ssn_itr(list_n)         =BB_lbfgsssn_itr.cons_ssn;
                A_ada_f1_ssn_itr(list_n)      =BB_lbfgsssn_itr.NewtonStep_BreakbyUnacceptedReduce;
                A_ada_f2_ssn_itr(list_n)      =BB_lbfgsssn_itr.NewtonStep_BreakbyDiffActiveSet;
                A_ada_time_used(list_n)       =BB_lbfgsssn_time_used;
                A_ada_normG(list_n)           =nrmG1;
                Diftheta0(list_n)             =diftheta;  Diftheta(list_n)=theta-theta1;  Diftheta1(list_n)=theta0-theta;
                Difnrm0(list_n)               =difnrm;    Difnrm(list_n)=nrmG-nrmG1   ;   Difnrm1(list_n)=nrmG0-nrmG;
                C_itr(list_n)                 =CADM_itr;
                C_time_used(list_n)           =CADM_time_used;
                P_itr(list_n)                 =pbfgs_itr;
                P_time_used(list_n)           =pbfgs_time_used;
                P_normG(list_n)               =nrmG2;
                
                Condnum(list_n)               =condnum;
                N(list_n) =matrix_n(mi); type(list_n)=casetype;
                le(list_n) =constr_n(1); ll(list_n) =constr_n(2); lu(list_n)=constr_n(3);
                
                CaseID(list_n)=string([num2str(casetype),'_',num2str(matrix_n(mi)),...
                    '_',num2str(constr_n(1)),'_',num2str(constr_n(2)),'_',num2str(constr_n(3))]);
                
%                     save("TestResult\temp"+experiments_date+".mat", ...
%                         "A_itr",    "A_time_used",      "A_normG",...
%                         "A_ada_itr","A_ada_time_used",  "A_ada_normG",...
%                         "B_itr",    "B_time_used",      "B_normG",...
%                         "D_itr",    "D_time_used",           ...
%                         "E_itr",    "E_time_used",      "E_normG",...
%                         "Condnum",...
%                         "Diftheta0", "Difnrm0",   "Diftheta","Difnrm",    "Diftheta1","Difnrm1" )
%                     load("TestResult\temp"+experiments_date+".mat");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
        end
    end
end
if test_N==1 %每一轮存一次结果。
    Test_Information=table(CaseID',type',N',le',ll',lu',...
        Newton_itr',   Newton_time_used',   Newton_normG',...
        A_Lbfgs_itr',  A_Lbfgs_time_used',  A_Lbfgs_normG',...
        A_ada_itr',    A_ada_time_used',    A_ada_normG', A_ada_ssn_itr', A_ada_f1_ssn_itr', A_ada_f2_ssn_itr',...
        C_itr',        C_time_used',           ...
        P_itr',        P_time_used',      P_normG',...
        Condnum',...
        Diftheta0', Difnrm0',   Diftheta',Difnrm',    Diftheta1',Difnrm1'  ,...
        'VariableNames',{'Case_Type_ID','type','N','le','ll','lu',...
        'SNM_Iterations',      'SNM_Time_Used',      'SNM_normG',...
        'LBFGS_Iterations',    'LBFGS_Time_Used',    'LBFGS_normG',...
        'LBFGSSSN_Iterations', 'LBFGSSSN_Iime_Used', 'LBFGSSSN_normG', 'LBFGSSSN_Newton_Iter','uR','uA',...
        'CADM_Iterations',     'CADM_Iime_Used',...
        'Pbfgs_Iterations',    'Pbfgs_Time_Used', 'Pbfgs_normG',...
        'CONDNUM',...
        'Dual_Fun_Diff_SNM_LbfgsSsn','Norm_Gred_Dif_SNM_LbfgsSsn',...
        'Dual_Fun_Diff_Lbfgs_LbfgsSsn','Norm_Gred_Dif_Lbfgs_LbfgsSsn',...
        'Dual_Fun_Diff_SNM_Lbfgs','Norm_Gred_Dif_SNM_Lbfgs'});
    path=[result_dir,'\SNM$BBLbfgs$BBLbfgsSsn$CADM$CONDNUM0','_ALLA',char(experiments_date),'.xls'];
    writetable(Test_Information, path, 'Sheet',datestr(now,30));    
end
fprintf('\n=============*==========*===========*=============*==========*==========')
fprintf('\n                        finish all the instances')
fprintf('\n=============*==========*===========*=============*==========*==========\n')

if exist('diary_on','var')==1, diary off;end %diary off




