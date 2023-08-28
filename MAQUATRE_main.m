% Conference Paper(ICGEC):A New Membrane Algorithm based on Quasi-Affine TRansformation Evolution for Optimization Problems
% Author: Wang Xiaopeng  Email: wangxp1993@163.com

clear;
clc;
format long;
fun_nums=12;
D=10;
Xmin=-100;
Xmax=100;
ps = 100;
iter_max=500;
runs=21;
%xbest 最优位置矩阵
xbest(runs,D) = inf;
xbest(:) = inf;
%fbest 最优值矩阵
fbest(fun_nums,runs) = inf;
fbest(:) = inf;
%error 矩阵
f_error(fun_nums,runs) = inf;
f_error(:) = inf;
%f_mean 均值矩阵
f_mean(fun_nums) = inf;
f_mean(:) = inf;
%f_median 中值矩阵
f_median(fun_nums) = inf;
f_median(:) = inf;
%f_std 标准差
f_std(fun_nums) = inf;
f_std(:) = inf;
Time = zeros(fun_nums,runs);

targetbest = [300;400;600;800;900;1800;2000;2200;2300;2400;2600;2700];

fhd=str2func('cec22_test_func');
fname = ['MAQUATRE_',num2str(D),'D.txt'];
f_out = fopen(fname,'wt');
fname_b_m_std_time = ['record_MAQUATRE_b_m_std_time',num2str(D),'D.txt'];
f_out_b_m_std_time = fopen(fname_b_m_std_time,'wt');
ftime = ['MAQUATRE_Time_',num2str(D),'D.txt'];
f_tout = fopen(ftime,'wt');

for i=1:fun_nums
    fun_num=i;
    disp(['fid:',num2str(fun_num)]);
    fprintf(f_out,'fid:%d\n',fun_num);
    fprintf(f_tout,'fid:%d\n',fun_num);
    for j=1:runs
        [gbest,gbestval,RecordT]= MAQUATRE_func(fhd,D,ps,iter_max,Xmin,Xmax,fun_num,j);
        xbest(j,:)=gbest;   %D维行向量
        fbest(i,j)=gbestval;%单一数值
        disp(['x[',num2str(gbest),']=',num2str(gbestval-targetbest(i),15)]);
        fprintf(f_out,'x%s\t[%s]=%s\n',num2str(j),num2str(gbest),num2str(gbestval-targetbest(i)));
        Time(i,j) = RecordT;
    end
    [bestval,best] = min(fbest(i,:));
    loc = xbest(best,:);
    disp(['Best[',num2str(loc),']=',num2str(bestval-targetbest(i),15)]);
    fprintf(f_out,'Best[%s]=%s\n',num2str(loc),num2str(bestval-targetbest(i)));

    f_error(i,:) = fbest(i,:) - targetbest(i);
    f_mean(i)=mean(f_error(i,:));
    f_median(i) = median(f_error(i,:));
    f_std(i) = std(f_error(i,:));
    disp(['mean[',num2str(i),']=',num2str(f_mean(i),15)]);
    fprintf(f_out,'mean[%s]=%s\n',num2str(i),num2str(f_mean(i)));
    disp(['median[',num2str(i),']=',num2str(f_median(i),15)]);
    fprintf(f_out,'median[%s]=%s\n',num2str(i),num2str(f_median(i)));
    disp(['std[',num2str(i),']=',num2str(f_std(i),15)]);
    fprintf(f_out,'std[%s]=%s\n',num2str(i),num2str(f_std(i)));
     
    for j = 1: runs
        disp(['Time[',num2str(j),']=',num2str(Time(i,j),15)]);
        fprintf(f_tout,'Time[%s]=\t%.15f\n',num2str(j),Time(i,j));
    end
    MeanT=mean(Time(i,:));
    disp(['meanT[',num2str(i),']=',num2str(MeanT,15)]);
    fprintf(f_tout,'MeanTime[%s]=\t%.15f\n',num2str(i),MeanT);
    
    best_mean_std(i,1)=i; %存到变量中
    best_mean_std(i,2)=bestval-targetbest(i);
    best_mean_std(i,3)=f_mean(i);
    best_mean_std(i,4)=f_std(i);
    best_mean_std(i,5)=MeanT;
    best_mean_std(i,6)=best;
    fprintf(f_out_b_m_std_time,'%s\n',num2str(best_mean_std(i,:)));

end
fclose(f_out);
fclose(f_tout);
fclose(f_out_b_m_std_time);
clear all;