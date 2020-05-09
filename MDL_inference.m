clear all

g = 10;
Tsamples = 100;
n = Tsamples - 1;

genes = [];

for i = 1:g
    genes = [genes,strcat('Gene',string(i))];
end
genes = cellstr(genes);

act_pred = [2,3,10;6,8,9;1,9,10;2,5,6;1,4,6;1,6,10;6,7,8;1,4,7;2,6,9;4,5,7];


num_datasets = 1000; 
K = 3;
tot = 0; for k = 1:K; tot=tot+size(combnk(1:g,k),1);end
pred = {};
for gene_id = 7
    
    genes(gene_id)
    disc_sample=0;
    
    L_M_combined = zeros(num_datasets,tot);
    L_N_combined= zeros(num_datasets,tot);
    for sample_id = 1:num_datasets
        run(strcat("100_samples_1000_datasets/random_network4_",string(sample_id)))
        data = [Gene1;Gene2;Gene3;Gene4;Gene5;Gene6;Gene7;Gene8;Gene9;Gene10];
        x = data(:,1:Tsamples-1);
        y = data(gene_id,2:Tsamples);
        if (sum(y) <= 1 || sum(y) >= n-1)
            %disp('all values are 0 or 1');
            disc_sample = disc_sample+1;
        end
        h = sum(y)/n;
        L0_N = n*(-h*log2(h)-(1-h)*log2(1-h)) + 0.5;
        L0_M = log2(Cml_calc(n))+ (1/2)*log2(pi/2) + log2(1 + log(g));

        sample_L_M = [];
        sample_L_N = [];

        for k = 1:K
            L_lambda = min(g, log2(nchoosek(g,k)) + log2(k+1) + log2(1+log(g)) );

            if ( L_lambda < L0_N + L0_M )
                H = combnk(1:g,k);
                L_M = zeros(1,size(H,1));
                L_N = zeros(1,size(H,1));
                for i = 1:size(H,1)
                    Xi = x(H(i,:),:);
                    [uniquerow, ~, rowidx] = unique(Xi.', 'rows');
                    noccurrences = accumarray(rowidx, 1);
                    w = length(noccurrences);
                    d = w;
                    ml = zeros(1,w);
                    ml_1 = zeros(1,w);
                    fval = zeros(1,w);
                    C_ml = zeros(1,w);
                    P_ml = zeros(1,w);
                    P_ml_log = zeros(1,w);
                    for l_idx = 1:w %2^k
                        ml(l_idx) = noccurrences(l_idx);
                        ml_1(l_idx) = sum(y(rowidx==l_idx));
            
                        var1 = y(rowidx==l_idx);
                        var2 = de2bi(ceil(ml(l_idx)/2));
           
                        if (ml_1(l_idx) < ml(l_idx) - ml_1(l_idx))
                            fval(l_idx) = 0;
                        elseif (ml_1(l_idx) > ml(l_idx) - ml_1(l_idx))
                            fval(l_idx) = 1;
                        elseif ((ml_1(l_idx) == ml(l_idx) - ml_1(l_idx)) && ( sum(var1.*10.^(size(var1, 2)-1:-1:0), 2) < sum(var2.*10.^(size(var2, 2)-1:-1:0), 2) ))
                            fval(l_idx) = 0;
                        else
                            fval(l_idx) = 1;
                        end
                        C_ml(l_idx) = Cml_calc(ml(l_idx));
                        %prop = sym(ml_1(l_idx)/ml(l_idx));
                        P_ml(l_idx) = ((ml_1(l_idx)/ml(l_idx))^(ml_1(l_idx)))*((1-ml_1(l_idx)/ml(l_idx))^(ml(l_idx)-ml_1(l_idx)));
                        %P_ml(l_idx) = vpa(prop^(ml_1(l_idx))*((1-prop)^(ml(l_idx)-ml_1(l_idx))),50);
                        P_ml_log(l_idx) = log2((ml_1(l_idx)/ml(l_idx)))*(ml_1(l_idx))+log2((1-ml_1(l_idx)/ml(l_idx)))*(ml(l_idx)-ml_1(l_idx));
                    end
                    P_y = prod(P_ml);
                    %P_y = vpa(prod(sym(P_ml)),50);
                    if(P_y <= 1e-20 ) 
                        P_y_log = sum(P_ml_log(P_ml~=1));
                        L_N(i) = -P_y_log + d/2;
                    else
                        L_N(i) = -log2(P_y) + d/2;
                    end
                    L_M(i) = sum(log2(C_ml))+ (w/2)*log2(pi/2) + L_lambda;
                end
                [val,id] = min(L_M+L_N);
                L0_M = L_M(id);
                L0_N = L_N(id);
                if(val==0)
                    disp('val is zero')
                end
            
                sample_L_M = [sample_L_M, L_M];
                sample_L_N = [sample_L_N, L_N]; 
            else  
                break;
            end
       
        end
        if (~isempty(sample_L_M) && ~isempty(sample_L_N))
            L_M_combined(sample_id,:) = sample_L_M;
            L_N_combined(sample_id,:) = sample_L_N;
        end
    end
    my_ind = find(~all(L_M_combined == 0,2));
    temp = sum(L_N_combined(my_ind,:),1)+mean(L_M_combined(my_ind,:),1);
    [sorted_val, sorted_id] = sort(temp);

    disp(strcat('Constant samples: ',string(disc_sample)))

    min_id = sorted_id(1);
    k = 1;
    H = combnk(1:g,k);


    while (size(H,1) < min_id)
        min_id = min_id - size(H,1);
        k = k+1; 
        H = combnk(1:g,k);    
    end

    H(min_id,:)
    pred(end+1) = {genes(H(min_id,:).')};
    disp('Rank of actual predictor set:')
    H = combnk(1:g,3);
    find(sorted_id-55 == find(ismember(H,act_pred(g,:),'rows')))
    disp('Most probable 3 gene predictor:')
    H(sorted_id(find(sorted_id>55, 1, 'first'))-55,:)
end
