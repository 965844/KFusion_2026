function nodes = adaptive_comm_update1_1(nodes,  i, iter, ub, lb, minKfusion, maxKfusion,mincons,Wt)
% ===============================================================
% 自适应通信更新函数 adaptive_comm_update
% ---------------------------------------------------------------
    maxf=2;
    dim=size(nodes(1).bestPos,2);

    
   % ---------- 2. 通信触发判断 ----------
    if iter == 1
    
        % ---------- 1. 构造双随机矩阵 ---------
        
        % 第一次迭代：强制通信
        [nodes(i).neighborBest,nodes(i).count] = get_neighbor_best(nodes, i,Wt,dim);  % 传节点 id
        nodes(i).kFusion = min(nodes(i).kFusion + 1, maxKfusion);
        nodes(i).didConsensus = true;
        nodes(i).consensusMetric = compute_consensus_metric(nodes, i, nodes(i).knownNeighbors); % 传邻居列表
    else
        currentMetric = compute_consensus_metric(nodes, i, nodes(i).knownNeighbors);
        % fprintf("邻居之间差距为%f,上一次邻居差距为%f\n",currentMetric,nodes(i).consensusMetric);
    
        if currentMetric <mincons
            nodes(i).inteval=nodes(i).inteval+1;
            nodes(i).didConsensus = false;
            nodes(i).kFusion = max(nodes(i).kFusion - 1, minKfusion);
        else
            % 共识恶化 -> 触发通信
            nodes(i).inteval=0;
            nodes(i).didConsensus = true;
            nodes(i).kFusion = min(nodes(i).kFusion + 1, maxKfusion);
            
        end
        nodes(i).consensusMetric = currentMetric;
    end
 
    if nodes(i).inteval >=nodes(i).maxinteval
        nodes(i).didConsensus=true;
    end
    % ---------- 3. 通信与共识更新 ----------
    if ~nodes(i).didConsensus
        % 无需通信，直接返回
        return;
    else
        % 获取最优解
        Wt=Wt^nodes(i).kFusion;
        [nodes(i).neighborBest,nodes(i).count] = get_neighbor_best(nodes, i,Wt,dim);  % 保证这里传 node id（i）
        % 通信 1：邻居最优解共识
        nodes(i).neighborBestMat = apply_Conse( ...
        Wt, nodes(i).neighborBest,ub,lb);
        % 计算共识
        
        nodes(i).consens=nodes(i).neighborBestMat(i,:);
        
    end

    nodes(i).kFusionHistory(iter) = nodes(i).kFusion;

end
function [neighborBest,count] = get_neighbor_best(nodes,i,Wt,dim)
% 提取邻居节点的最优解集合
   neighborBest = [];
   count=0;
   nodes(i).Neig=sum(Wt(i,:)~=0);
   for idx = 1:size(Wt,2)
       if Wt(i,idx) ~=0
           count=count+1;
           neighborBest=[neighborBest;nodes(idx).bestPos];
       else
           neighborBest=[neighborBest;nodes(i).bestPos];
       end
     
   end
end
function metric = compute_consensus_metric(nodes, i, neighbors)
% 计算当前节点与邻居最优解之间的差异指标
    if isempty(neighbors)
        metric = 0;
        return;
    end
    selfSol = nodes(i).consens;
    diffs = zeros(length(neighbors), 1);
    for k = 1:length(neighbors)
        diffVal = norm(selfSol - nodes(neighbors(k)).consens);
        diffs(k) = diffVal;
    end
    metric = mean(diffs);
end

function resultMat = apply_Conse( W, M,ub,lb)
    resultMat= max(min(W * M, ub), lb);
end
