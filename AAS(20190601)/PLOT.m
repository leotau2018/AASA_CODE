%% Draw residual(e) and time(t) on E5 and E6 with n=1000 or n=5000, nr=500
%% And Real Case E7 and E8
%--------------------------------------------------------------------------


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% E5 or E6
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% select figure parameters
% input data
clear
plot_example_type = 6;
plot_example_dim  = 2000;
table_path = '.\PLOT_data\';

% transport parameters
if plot_example_dim == 1000
    dime = '1000';
else
    dime = '2000';
end
table_name = [table_path,'diary_',dime,'_E5E6.xlsx'];
if plot_example_type == 5
    E_type='E5'; E_typeinarticle = 'E3';
else
    E_type='E6'; E_typeinarticle = 'E4';
end

%% read data
data = xlsread(table_name,[E_type,dime,'_ISNM']);
e_ISNM = log10(data(:,3));
t_ISNM = data(:,4);

data = xlsread(table_name,[E_type,dime,'_AASA(LBFGS)']);
e_LBFGS     = log10(data(:,2));
t_LBFGS     = data(:,4);

data = xlsread(table_name,[E_type,dime,'_AASA(Adaptive)']);
e_SNN_LBFGS = log10(data(:,2));
t_SNN_LBFGS = data(:,4);
flag_SNN_LBFGS = data(:,5);


%% plot figures
%% -----------------------------------------------
% plot e 
figure
% subplot(1,2,1)
plot(e_SNN_LBFGS,'m-s','MarkerIndices',[1:10:length(e_SNN_LBFGS)-5,length(e_SNN_LBFGS)-4:length(e_SNN_LBFGS)]);hold on;
plot(e_LBFGS,'c:o','MarkerIndices',[1:30:length(e_LBFGS)-5,length(e_LBFGS)-4:length(e_LBFGS)])
plot(e_ISNM,'r-.d','MarkerIndices',[1:40:length(e_ISNM)-5,length(e_ISNM)-4:length(e_ISNM)])
% y_tick = {'10^{-8}','10^{-6}','10^{-4}','10^{-2}','10^0','10^2','10^{4}'};
y_tick = change_y_tick;%( gca);
set(gca, 'YtickLabel',y_tick);
xlim([1,max([length(e_SNN_LBFGS),length(e_LBFGS),length(e_ISNM)])+100]);
grid on;
time_SNN_LBFGS = ['AASA(Adaptive)(',num2str(max(t_SNN_LBFGS)),'s)'];
time_LBFGS     = ['AASA(L-BFGS)(',num2str(max(t_LBFGS)),'s)'];
time_ISNM      = ['AASA(L-BFGS)(',num2str(max(t_ISNM)),'s)'];
legend(time_SNN_LBFGS, time_LBFGS, time_ISNM)
title([E_typeinarticle,': n=',dime,' n_r=500'])
xlabel('Iteration number')
ylabel('Residuals')

% emphasis newton step 
% plot(e_SNN_LBFGS,'s','MarkerIndices',find(flag_SNN_LBFGS>1), 'MarkerFaceColor','m','MarkerEdgeColor','m')
% legend(time_SNN_LBFGS, time_LBFGS, time_ISNM, 'Newton Step')
%% -----------------------------------------------
%%%% plot t
figure
% subplot(1,2,2)
plot(t_SNN_LBFGS,'m-s','MarkerIndices',[1:10:length(t_SNN_LBFGS),length(t_SNN_LBFGS)]);hold on;
plot(t_LBFGS,'c:o',    'MarkerIndices',[1:30:length(t_LBFGS),length(t_LBFGS)])
plot(t_ISNM,'r-.d',    'MarkerIndices',[1:40:length(t_ISNM),length(t_ISNM)])
xlim([1,max([length(e_SNN_LBFGS),length(e_LBFGS),length(e_ISNM)])+100]);
grid on;
itr_SNN_LBFGS = ['AASA(Adaptive)(',num2str(length(t_SNN_LBFGS)),')'];
itr_LBFGS     = ['AASA(L-BFGS)(',num2str(length(t_LBFGS)),')'];
itr_ISNM      = ['AASA(L-BFGS)(',num2str(length(t_ISNM)),')'];
legend(itr_SNN_LBFGS, itr_LBFGS, itr_ISNM, 'Location','best')
legend('AASA(Adaptive)', 'AASA(L-BFGS)', 'AASA(L-BFGS)', 'Location','best')
title([E_typeinarticle,': n=',dime,' n_r=500'])
xlabel('Iteration number')
ylabel('CPU time (s)')




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Draw residual(e) on Real Case E7 and E8
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% input data
clear
plot_example_type = 6;
plot_example_dim  = 2000;
table_path = '.\PLOT_data\';
table_name = [table_path,'diary_RE7RE8.xlsx'];
%% read data
n=792; %792£¬   1187
data = xlsread(table_name,['CADM',num2str(n)]);
e_CADM = log10(data(1:end-1,2));
t_CADM = data(end,2);

data = xlsread(table_name,['Pbfgs',num2str(n)]);
e_PBFGS = log10(data(1:end-1,2));
t_PBFGS = data(end,2);

data = xlsread(table_name,['ISNM',num2str(n)]);
e_ISNM = log10(data(1:end-1,2));
t_ISNM = data(end,2);

data = xlsread(table_name,['AASA(LBFGS)',num2str(n)]);
e_LBFGS     = log10(data(1:end-1,2));
t_LBFGS     = data(end,2);

data = xlsread(table_name,['AASA(Adaptive)',num2str(n)]);
e_SNN_LBFGS = log10(data(1:end-1,2));
t_SNN_LBFGS = data(end,2);
% flag_SNN_LBFGS = data(:,5);

%% -----------------------------------------------
%%%% plot t
figure
plot(e_CADM,'b--.','MarkerIndices',[1:30:length(e_CADM),length(e_CADM)],'MarkerSize',15);hold on; %semilogy(e_CADM)
% ylim([-8,4]);
plot(e_SNN_LBFGS,'m-s');%semilogy(e_SNN_LBFGS)
plot(e_LBFGS,'c:o','MarkerIndices',[1:10:length(e_LBFGS),length(e_LBFGS)])
plot(e_ISNM,'r-.d')
plot(e_PBFGS,'g:p','MarkerIndices',[1:10:length(e_PBFGS),length(e_PBFGS)])
% y_tick = {'10^{-8}','10^{-6}','10^{-4}','10^{-2}','10^0','10^2','10^{4}'};
y_tick = change_y_tick;
set(gca, 'YtickLabel',y_tick);
grid on;
annSNN_BFGS=annotation('textarrow',[0.35 0.15],[0.22 0.22],...
     'String',['CPU:',num2str(t_SNN_LBFGS),'(s)'],'Color','m');annSNN_BFGS.FontSize = 12;
annLBFGS=annotation('textarrow',   [0.35 0.18],[0.34 0.34],...
     'String',['CPU:',num2str(t_LBFGS),'(s)'],'Color','c');annLBFGS.FontSize = 12;
annISNM=annotation('textarrow',    [0.35 0.15],[0.16 0.16],...
     'String',['CPU:',num2str(t_ISNM),'(s)'],'Color','r');annISNM.FontSize = 12;
annCADM=annotation('textarrow',[0.77 0.77],[0.35 0.27],...
     'String',['CPU:',num2str(t_CADM),'(s)'],'Color','b');annCADM.FontSize = 12;
annPBFGS=annotation('textarrow',[0.35 0.28],[0.28 0.28],...
     'String',['CPU:',num2str(t_PBFGS),'(s)'],'Color','g');annPBFGS.FontSize = 12;
% text(length(e_CADM)-200, e_CADM(end)-0.5, ' 1042.3(s) \uparrow ' ,  'Color','b','FontSize',12);
% text(length(e_SNN_LBFGS),e_SNN_LBFGS(end),' \leftarrow 101.1(s)' ,  'Color','m','FontSize',12);
% text(length(e_LBFGS),    e_LBFGS(end),    ' \leftarrow 152.7(s)' ,  'Color','c','FontSize',12);
% text(length(e_ISNM)-30,  e_ISNM(end)+0.5, ' \downarrow 354.7(s)' ,  'Color','r','FontSize',12);
% text(length(e_PBFGS),    e_PBFGS(end),    ' \leftarrow 828.4(s) ',  'Color','g','FontSize',12);
legend('ADM','AASA(Adaptive)','AASA(L-BFGS)','ISNM','P-BFGS');
title('n=792, n_r=n')
if n==1187, title('n=1187, n_r=400');end
xlabel('Iteration number')
ylabel('Residuals')



%% change YTick %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------
function new_y_tick = change_y_tick
mygca = gca;
new_y_tick = {};
for i = 1:length(mygca.YTick)
    new_y_tick{i} = ['10^{',num2str(mygca.YTick(i)),'}'];
end
end