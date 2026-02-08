function nodes = node_evolution_update1(nodes, i,   ub, lb,  R, A,xopt)
    node_i = nodes(i);
    pop = node_i.population;
    fit = node_i.fitness;
    popSize = size(pop,1);
    l=nodes(i).bestPos;


    F_base=0.9;%基础变异强度

    % -------- 计算群体多样性 --------
    center = mean(pop,1);
    diversity = mean(vecnorm(pop - center, 2, 2));
    
    % 根据种群多样性对强度进行调整
    if diversity < 3
        % 多样性不足 → 加强探索
        F = min(1.3, F_base * 1.5);
    else
        % 多样性充足 → 稳定搜索··
        F = F_base;
    end
    
    %% 第一种进化方式
    for j = 1:popSize

        r = randsample(popSize, 3);
        v = pop(r(1),:) + F * (pop(r(3),:) - pop(r(2),:));
        v = max(min(v, ub), lb);
        fv = benchmark_functions(nodes(i), v,  R, A,xopt);

        if fv < fit(j)
            pop(j,:) = v;
            fit(j) = fv;
        end
    end
    
 

    nodes(i).population = pop;
    nodes(i).fitness = fit;


    %% 找出更新后的最优解和值
    [bf, idx2] = min(nodes(i).fitness);
    nodes(i).bestFitness = bf;
    nodes(i).bestPos = nodes(i).population(idx2,:);

    dis=norm(nodes(i).bestPos - l);
    if dis==0
        nodes(i).bestcount=nodes(i).bestcount+1;
    else
        nodes(i).bestcount=max(0,nodes(i).bestcount-1);
    end
    
    % fprintf("策略2结点%d本次更新最优解差距为%d\n",i,dis);
    
    
    


end
