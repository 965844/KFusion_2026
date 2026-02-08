function nodes = adaptive_comm_update1_2(nodes,  i, iter, ub, lb, minKfusion, maxKfusion,mincons,Wt)
% ===============================================================
% 自适应通信更新函数 adaptive_comm_update
% ---------------------------------------------------------------
    dim=100;
    
 
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
        % 计算共识
        nodes(i).consens=max(min(Wt(i,:)*nodes(i).neighborBest, ub), lb);
    end
    nodes(i).kFusionHistory(iter) = nodes(i).kFusion;
end
function [neighborBest,count] = get_neighbor_best(nodes,i,Wt,dim)
% 提取邻居节点的最优解集合
   neighborBest = [];
   count=0;
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

