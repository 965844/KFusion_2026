clc; clear;

%% ---------- 基本参数 ----------
numNodes = 20;

funid=10;                                                       

popSize  = 300;
dim      =100;
maxIter  =2000;
lb = -100;  
ub=100;
c_cons   =0.2;
injectRate = 0.1;

minconus=0.5;




%1随机图 2环图 3全连接图   邻接矩阵


% 用于判定通信次数
count=zeros(1,maxIter+1);

% 根据函数id获取数据

xopt=getX(funid);
R=getR(funid);
W=getW(funid);
A=getA(funid);


con_change=[];


glo_vla=zeros(1,maxIter);
cons_vla=zeros(1,maxIter);

minKfusion = 0;
maxKfusion = 20; 

function W=getW(id)
    strpath="data\W_F"+string(id);
    W=load(strpath);
end
function X=getX(id)
    strpath="data\xopt_F"+string(id);
    X=load(strpath);
    
end
function R=getR(id)
    strpath="data\R_F"+string(id);
    R=load(strpath);
end
function A=getA(id)
    strpath="data\A_F"+string(id);
    A=load(strpath);
end




%% ---------- 初始化各节点 ----------
nodes = struct([]);
for i = 1:numNodes
    
    pop = lb + (ub-lb) .* rand(popSize, dim);
    nodes(i).funid=funid;

   
    
    nodes(i).count=0;
    nodes(i).id = i;
    fit=benchmark_functions(nodes(i),pop,R,A,xopt);
    nodes(i).population = pop;
    nodes(i).fitness    = fit;

    [bestFit,bestIdx] = min(fit);
    nodes(i).bestPos     = pop(bestIdx,:);
    nodes(i).bestFitness = bestFit;
    nodes(i).bestFitnessHistory = bestFit;

    nodes(i).hisFit=bestFit;
    nodes(i).hisPos= pop(bestIdx,:);

    nodes(i).inteval=0;             %连续多少次没通信
    nodes(i).maxinteval=2;
    
    nodes(i).bestPosHistory=[];     % 每次迭代最优解
    nodes(i).kFusion     = 1;       % 初始通信强度
    nodes(i).W           = [];      % 暂无权重矩阵（后续学习）
    nodes(i).neighborBest = zeros(numNodes,dim);     % 尚无邻居信息
    nodes(i).consensusMetric = Inf;
    nodes(i).didConsensus = false;  %是否通信
    nodes(i).neighborBestMat = [];
    nodes(i).U=[];
    nodes(i).consens=nodes(i).bestPos;

    nodes(i).bestcount=0;

    nodes(i).sigma_g=0.5;
  
    
    nodes(i).kFusionHistory= zeros(maxIter,1);
    
    directNbrs = find(W(i,:) ~= 0);
    expandedList = unique([i, directNbrs]);

    nodes(i).knownNeighbors = expandedList;%存放邻居的编号
    nodes(i).numNeighbors=length(expandedList);

    neighbors = find(W(i,:) ~= 0);   % 从邻接矩阵提取邻居
    numNeighbors = length(neighbors)+1;

    nodes(i).evolutionid=1;

    nodes(i).jiange=[5,10];
    nodes(i).celue=10;
    nodes(i).Neig=0;
    
end
nnn=0;
count1=0;

c = parcluster('local');
c.NumWorkers = 20;
saveProfile(c);

if isempty(gcp('nocreate'))
    parpool(c.NumWorkers);
end
kf_max=zeros(1,maxIter);
for iter = 1:maxIter
   
    nodes_new = nodes;
    sum_global = 0;
    kfu=0;
    kFusion_local = zeros(numNodes,1);
    nodeMean = zeros(numNodes, dim);
    nodes_tmp = nodes;
    %多线程跑
    parfor kk=1:numNodes
      
      nodes_tmp=nodes;
       %% 
       % 间隔通信
       nodes_tmp = adaptive_comm_update1_1(nodes_tmp, kk, iter, ub,lb, minKfusion, maxKfusion,minconus,W);
       %%
       %%自适应进化策略选择
       if  nodes_tmp(kk).bestcount<=1
                nodes_tmp(kk).celue=1;
                nodes_tmp = node_evolution_update(nodes_tmp, kk, iter,  c_cons, ub, lb,maxIter,R,A,xopt);
       elseif nodes_tmp(kk).bestcount>1 && nodes(kk).bestcount<=5
                nodes_tmp(kk).celue=2;
                nodes_tmp = node_evolution_update1(nodes_tmp, kk,  ub, lb,R,A,xopt);
       elseif nodes_tmp(kk).bestcount>5
                nodes_tmp(kk).celue=3;
                nodes_tmp = node_evolution_update2(nodes_tmp, kk, iter,maxIter,0.3, ub, lb,R,A,xopt);
       end
       
        %%
        %更新本地最优
            [bf, idx2] = min(nodes_tmp(kk).fitness);
            nodes_tmp(kk).bestFitness = bf;
            nodes_tmp(kk).bestPos = nodes_tmp(kk).population(idx2,:);
        %% 注入：强化局部收敛
            [~, idxDesc] = sort(nodes_tmp(kk).fitness, 'descend');
            kInject = max(1, round(injectRate * popSize));
            worstIdx = idxDesc(1:kInject);
        %%
        %开始注入
            part1 = nodes_tmp(kk).bestPos + 0.05 * randn(kInject, dim);
            part2 = nodes_tmp(kk).consens + 0.05 * randn(kInject, dim);
            nodes_tmp(kk).population(worstIdx,:) = 0.5 * part1 + 0.5 * part2;
    
            nodes_tmp(kk).population(worstIdx,:) = max(min(nodes_tmp(kk).population(worstIdx,:), ub), lb);
    
            nodes_tmp(kk).fitness(worstIdx) = ...
                arrayfun(@(r) benchmark_functions(nodes_tmp(kk), nodes_tmp(kk).population(r,:),  R, A,xopt), worstIdx);
        %%
        %%计算共识差距和评价
        nodeMean(kk,:) = nodes_tmp(kk).consens;
        % 用 bestFitness 统计全局值
        sum_global = sum_global + benchmark_functions(nodes_tmp(kk),nodes_tmp(kk).consens,R,A,xopt);
        nodes_new(kk) = nodes_tmp(kk);
    
    end
        
        nodes = nodes_new;
        kf_max(iter)=mean([nodes.kFusion]);
        glo_vla(iter) = sum_global / numNodes;
        % 共识距离
        meanCons = mean(nodeMean, 1);
        cons_vla(iter) = sum( sqrt(sum((nodeMean - meanCons).^2, 2)) );
        fprintf("第%d次迭代最好的结果是%.4e,目前共识为%.4e,\n",iter,glo_vla(iter),cons_vla(iter));
   

end
fprintf("\n结果为%f 通信次数%d\n ",sum_global / numNodes,count(iter+1));




%探索邻居 变化显示
figure; 
plot(kf_max,'LineWidth',1.5);
xlabel('Iteration'); ylabel('value');
title('探索邻居(kFusion)变化');
grid on;


%全局目标值 变化显示
figure; 
plot(glo_vla,'LineWidth',1.5);
xlabel('Iteration'); ylabel('value');
title('全局目标值变化');
grid on;

%共识 变化显示
figure; 
plot(cons_vla,'LineWidth',1.5);
xlabel('Iteration'); ylabel('value');
title(numNodes,'结点共识变化');
grid on;

% %通信次数显示
% figure; 
% plot(count,'LineWidth',1.5);
% xlabel('Iteration'); ylabel('value');
% title(numNodes,'结点通信次数');
% grid on;


