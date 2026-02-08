function nodes = node_evolution_update(nodes, i, iter, c_cons, ub, lb, maxIter,   R, A,xopt)
    node_i = nodes(i);
    pop = node_i.population;
    fit = node_i.fitness;
    popSize = size(pop,1);
    l=nodes(i).bestPos;


  

    % 稳定动态 F
    F = 0.4- 0.3 * (iter / maxIter);
    F=F*2;

    alpha = c_cons * (iter / maxIter);
    alpha = min(alpha, 0.7);
    %% 第一种进化方式
    for j = 1:popSize

        r = randsample(popSize, 3);

        v = pop(r(1),:) + F * (nodes(i).consens - pop(r(2),:));
        % 凸组合共识
        new_x = (1-alpha) * v + alpha * nodes(i).consens;
        
        new_x = max(min(new_x, ub), lb);

        fnew = benchmark_functions(nodes(i), new_x,  R, A,xopt);

        if fnew < fit(j)
            pop(j,:) = new_x;
            fit(j) = fnew;
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
   
    % fprintf("策略1 结点%d本次更新最优解差距为%d\n",i,dis);

    


end
