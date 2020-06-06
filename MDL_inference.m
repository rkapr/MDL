 
clear

minver = '9.7';
if ~strncmp(version,minver,1)
  warning('MYTEST:VERCHK',...
      'Script verified on MATLAB 2019b v9.7, currently running v%s',minver)
end

format compact

g = 10; %Number of genes in the network
Tsamples = 10; %length of time series

num_datasets = 10; %Number of datasets, each dataset is a time series of 
                    %length Tsamples
n = num_datasets*(Tsamples - 1);
                    
% dfile =strcat('Output_Ts',string(Tsamples),'_Ns',string(num_datasets),'.txt');
% if exist(dfile, 'file') ; delete(dfile); end
% diary(dfile)
% 
% 
% diary on

load('all_data.mat')
gene_matrix2 = reshape(gene_matrix,10,100,1000);


gene_list = 1:10;

% Actual Predictors
act_pred = [2,3,10;6,8,9;1,9,10;2,5,6;1,4,6;...
    1,6,10;6,7,8;1,4,7;2,6,9;4,5,7];

K = 3; % Maximum number of allowed predictors per gene
H_numEl_arr=zeros(1,K); % Array containing number of possible predictor
                        % set for each k = 1:K
for k = 1:K; H_numEl_arr(k)=size(combnk(1:g,k),1);end
tot = sum(H_numEl_arr); % Total number of predictor sets
pred = {}; % Cell array to store predictor gene sets for each gene
fval_cell = {};
uniquerow_cell = {};
network_cost = zeros(g,tot);
% Array of gene names
genes = [];
for i = 1:g
    genes = [genes,strcat('Gene',string(i))];
end
genes = cellstr(genes);

% For SF calculations:
% maxval = zeros(g,num_datasets);
% SF_counts = zeros(g,tot);

for gene_id = gene_list
    
    disp(genes(gene_id))
    disc_sample=0; % Varible to store number of datasets with constant gene
                   % expression
    
    % Matrices to store model and noise codelength for each predictor set
    % and each sample
    L_M_combined = [];
    L_N_combined= [];
    
        y=gene_matrix2(gene_id,2:Tsamples,1:num_datasets);
        y=reshape(y,1,[]);
        x=gene_matrix2(:,1:Tsamples-1,1:num_datasets);
        x=reshape(x,10,[]);
        if (sum(y) < 1  || sum(y) > n-1)
            % Removing this improves detection for smaller # datasets
            %disp('all values are 0 or 1');
            disc_sample = 1;
            %continue;
        end
        h = sum(y)/n;
        L0_N = n*(-h*log2(h)-(1-h)*log2(1-h)) + 0.5;
        L0_M = log2(Cml_calc(n))+ (1/2)*log2(pi/2) + log2(1 + log(g));

        for k = 1:K
            L_lambda = min(g, log2(nchoosek(g,k)) + log2(k+1) + log2(1+log(g)) );

            if ( L_lambda < L0_N + L0_M )
                H = combnk(1:g,k);
                L_M = zeros(1,size(H,1));
                L_N = zeros(1,size(H,1));
                
                for i = 1:size(H,1)
                    Xi = x(H(i,:),:);
                    [uniquerow, ~, rowidx] = unique(Xi.', 'rows');
                    uniquerow_cell(gene_id,k,i)={uniquerow};
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
                        var2 = de2bi_code(ceil(ml(l_idx)/2));
           
                        if (ml_1(l_idx) < ml(l_idx) - ml_1(l_idx))
                            fval(l_idx) = 0;
                        elseif (ml_1(l_idx) > ml(l_idx) - ml_1(l_idx))
                            fval(l_idx) = 1;
                        elseif ((ml_1(l_idx) == ml(l_idx) - ml_1(l_idx)) && ( sum(var1.*10.^(size(var1, 2)-1:-1:0), 2) < sum(var2.*10.^(size(var2, 2)-1:-1:0), 2) ))
                            fval(l_idx) = 0;
                        else
                            fval(l_idx) = 1;
                        end
                        fval_cell(gene_id,k,i) = {fval};
                        C_ml(l_idx) = Cml_calc(ml(l_idx));
                        P_ml(l_idx) = ((ml_1(l_idx)/ml(l_idx))^(ml_1(l_idx)))*((1-ml_1(l_idx)/ml(l_idx))^(ml(l_idx)-ml_1(l_idx)));
                        P_ml_log(l_idx) = log2((ml_1(l_idx)/ml(l_idx)))*(ml_1(l_idx))+log2((1-ml_1(l_idx)/ml(l_idx)))*(ml(l_idx)-ml_1(l_idx));
                    end
                    P_y = prod(P_ml);
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
            
                L_M_combined = [L_M_combined, L_M];
                L_N_combined = [L_N_combined, L_N]; 
            end

       
        end
        % For SF calculations:
%         temp1 = [sample_L_M];
%         temp2 = [sample_L_N];
%         kk=convhull(temp1,temp2);
%         temp3=temp1(kk);
%         temp4=temp2(kk);
%         slope_val = (temp4(2:end) - temp4(1:end-1))./(temp3(2:end)-temp3(1:end-1));
%         if ~isempty(max(kk(slope_val<-1)))
%             %maxval(gene_id,sample_id) = max(kk(slope_val<-1));
%             [SF_sample_counts,SF_id]=groupcounts(kk(slope_val< -1));
%             SF_counts(gene_id,SF_id) = SF_counts(gene_id,SF_id)+SF_sample_counts.';
%         end
    my_ind = find(~all(L_M_combined == 0,2));
    if(~all(my_ind == 0,2))
        temp = sum(L_N_combined(my_ind,:),1)+mean(L_M_combined(my_ind,:),1);
        [sorted_val, sorted_id] = sort(temp);
        
        disp('Actual predictor set:')
        disp(act_pred(gene_id,:))
        
        min_id = sorted_id(1)-[0,cumsum(H_numEl_arr(1:end-1))];
        k = find(min_id>0,1,'last');
        H = combnk(1:g,k);
        disp('Estimated predictor set:')
        disp(H(min_id(k),:))
        pred(end+1) = {genes(H(min_id(k),:).')};
        disp('Rank of actual predictor set:')
        H = combnk(1:g,3);
        disp(find(sorted_id-55 == find(ismember(H,act_pred(gene_id,:),'rows'))))
        disp('Most probable 3 gene predictor set:')
        disp(H(sorted_id(find(sorted_id>55, 1, 'first'))-55,:))
        %cell2mat(fval_cell(gene_id,k,min_id(k)))
        %cell2mat(uniquerow_cell(gene_id,k,min_id(k)))
        
        % For SF calculations:
        %     peaks = find(SF_counts(gene_id,:)>0);
        %     min_id = peaks(end)-[0,cumsum(H_numEl_arr(1:end-1))];
        %     kval = find(min_id>0,1,'last');
        %     H = combnk(1:g,kval);
        %     disp('Estimated predictor set SF:')
        %     disp(H(min_id(kval),:))
        
        %     [SF_res_len,SF_res_ids] = sort(temp(peaks(end-1):peaks(end)));
        %     if (all(SF_res_ids+peaks(end-1)>55))
        %
        %         disp('Other likely 3 gene predictor sets SF (increasing codelengths):')
        %         disp(H(SF_res_ids+peaks(end-1)-55,:))
        %     end
        sorted_id_genes(gene_id,:) = sorted_id;
    else
        disp(strcat('All constant data: ',num2str(disc_sample)))
        disp('Estimated predictor set: NULL')
    end
    network_cost(gene_id,:) = temp;
    % run("est_PBN")
end
%diary off
% save('Ts100_Ns1000')
% save('Ts100_Ns1000/uniquerow_cell_Ts100_Ns1000','uniquerow_cell')
% save('Ts100_Ns1000/fval_cell_Ts100_Ns1000','fval_cell')
% save('Ts100_Ns1000/sorted_id_genes_Ts100_Ns1000','sorted_id_genes')
% save(strcat('Ts',string(Tsamples),'_Ns',string(num_datasets)))

valid_mat = zeros(g,tot);
cumsm_H_numEl = cumsum(H_numEl_arr);

for gene_id = 1:10
    if sorted_id_genes(gene_id,1)==0
        continue;
    end
valid_mat(gene_id,sorted_id_genes(gene_id,1))=1;
min_id = sorted_id_genes(gene_id,1)-[0,cumsm_H_numEl(1:end-1)];
k = find(min_id>0,1,'last');

H = combnk(1:g,k);
pred_id = H(min_id(k),:);
for k1=k+1:3
    H1 = combnk(1:g,k1);
    valid_mat(gene_id,cumsm_H_numEl(k1-1)+find(sum(ismember(H1,pred_id),2)==k)) = 1;
end
end


format loose
