function nodes = node_evolution_update2(nodes, i, iter, maxIter, c_cons, ub, lb, R, A, xopt)

    node_i  = nodes(i);
    pop     = node_i.population;
    fit     = node_i.fitness;
    popSize = size(pop,1);
    dim     = size(pop,2);

    lastBest = nodes(i).bestPos;

    % =========================
    % 1) 计算多样性 & 自适应步长 sigma
    % =========================
    center = mean(pop,1);
    div = mean(vecnorm(pop - center, 2, 2));  % 一个标量多样性指标

    % per-dim尺度：避免不同维度尺度差异导致某些维度“基本不动”
    sigma_vec = std(pop, 0, 1);
    sigma_vec = max(sigma_vec, 1e-12);

    % 全局步长（随迭代略收缩，但不会收缩到0）
    % 你可以把 0.2/0.02 调成适合你问题尺度的量级
    sigma_g = (0.20 * (1 - iter/maxIter) + 0.02);

    % 多样性太小 -> 适当放大步长，防止塌缩后“差分也很小”
    if div < 3
        sigma_g = sigma_g * 1.8;
    end

    % 共识融合强度
    alpha = min(0.5, c_cons * (iter/maxIter));  % 比策略1温和一些

    xbest = nodes(i).bestPos;
    xcons = nodes(i).consens;

    % =========================
    % 2) ES 采样：高斯局部 + 俄罗斯大跳（Cauchy）
    % =========================
    improved = 0;

    % 交叉概率：保留部分原基因，避免全维一起乱跳导致破坏性过大
    CR = 0.3 + 0.4*(1 - iter/maxIter);  % 前期更敢换，后期更保守

    for j = 1:popSize
        x = pop(j,:);

        % 以一定概率围绕 best 或 consensus 采样
        if rand < 0.6
            base = (1-alpha)*xbest + alpha*xcons;   % 强利用：靠近好解但仍带一点共识
        else
            base = x;                                % 保留个体自身探索
        end

        % --- 选择扰动分布 ---
        if rand < 0.8
            % 高斯扰动：精细搜索
            step = randn(1,dim) .* (sigma_g .* sigma_vec);
        else
            % 柯西扰动：偶尔俄罗斯大跳，专治“跳不出小坑”
            step = tan(pi*(rand(1,dim)-0.5)) .* (0.6*sigma_g .* sigma_vec);
        end

        y = base + step;

        % Binomial crossover：只替换部分维度
        mask = rand(1,dim) < CR;
        if ~any(mask), mask(randi(dim)) = true; end
        trial = x;
        trial(mask) = y(mask);

        % 边界处理（你现在是硬clip，这里保持一致）
        trial = max(min(trial, ub), lb);

        ftrial = benchmark_functions(nodes(i), trial, R, A, xopt);

        if ftrial < fit(j)
            pop(j,:) = trial;
            fit(j) = ftrial;
            improved = improved + 1;
        end
    end

    % =========================
    % 3) 1/5 success rule 调整“下一代步长状态”
    %    （用 bestcount/字段把 sigma 保存下来更好）
    % =========================
    successRate = improved / popSize;

    % 把 sigma 存到 nodes(i).sigma_g 里实现跨代记忆
    if ~isfield(nodes(i), 'sigma_g') || isempty(nodes(i).sigma_g)
        nodes(i).sigma_g = sigma_g;
    end

    if successRate > 0.2
        nodes(i).sigma_g = min(nodes(i).sigma_g * 1.15, 2.0);
    else
        nodes(i).sigma_g = max(nodes(i).sigma_g / 1.15, 1e-6);
    end

    % =========================
    % 4) 停滞触发：部分重启（替换最差 20%）
    % =========================
    nodes(i).population = pop;
    nodes(i).fitness = fit;

    [bf, idxBest] = min(fit);
    nodes(i).bestFitness = bf;
    nodes(i).bestPos = pop(idxBest,:);

    dis = norm(nodes(i).bestPos - lastBest);
    if dis == 0
        nodes(i).bestcount = nodes(i).bestcount + 1;
    else
        nodes(i).bestcount = 0;
    end

    % 连续不动很多次 -> 重启最差个体
    stallK = 25;  % 可以调：20~50
    if nodes(i).bestcount >= stallK
        % 找最差的 20%
        [~, ord] = sort(nodes(i).fitness, 'descend');
        k = max(1, round(0.2 * popSize));
        worstIdx = ord(1:k);

        % 以 best 为中心扩散 + 少量全局随机注入
        for t = 1:k
            if rand < 0.7
                % 围绕 best 做较大扰动（扩大探索半径）
                z = xbest + randn(1,dim) .* (3.0 * nodes(i).sigma_g .* sigma_vec);
            else
                % 少量全局随机（防止全陷在一个区域）
                z = lb + rand(1,dim) .* (ub - lb);
            end
            z = max(min(z, ub), lb);
            fz = benchmark_functions(nodes(i), z, R, A, xopt);

            pop(worstIdx(t),:) = z;
            fit(worstIdx(t)) = fz;
        end

        nodes(i).population = pop;
        nodes(i).fitness = fit;

        [bf, idxBest] = min(fit);
        nodes(i).bestFitness = bf;
        nodes(i).bestPos = pop(idxBest,:);

        nodes(i).bestcount = 0;  % 重启后清零
        
        
        
    end
    dis=norm(nodes(i).bestPos - xbest);
    if dis==0
        nodes(i).bestcount=nodes(i).bestcount+1;
    else
        nodes(i).bestcount=max(0,nodes(i).bestcount-0.2);
    end
       
    % fprintf("策略3 结点%d本次更新最优解差距为%d\n",i,dis);
   
end
