function W = getW(tt, p,maxInDegree, w,q,n,eta,tau, inter_type, monit_type,b,A,mu,Sigma,A_oc,mu_oc,Sigma_oc,doValue)
%UNTITLED calculate the monitoring statictics
% inpout: tt: conut index
%         p:number of nodes
%         maxInDegree: d
%         w,q eta: window parameter 
%         n: max run length
%         tau: time when change occurs 
%         inter_type, monit_type: methods
%         b: thresherhold. set to Inf
%         A mu Sigma: pre-change parameter for DAG
%         A_oc mu_oc Sigma_oc: post-change parameter for DAG
%         doValue: intervetnion value.
%output: W: monitoring statistics
seqN = zeros(1, w);
indices = randperm(w, ceil(q^(eta)));
seqN(indices) = 1;

seqN = repmat(seqN,1,ceil(n/w)); 
seqN = seqN(1:n);

%%%%%%%
doSet = -1*ones(1,n); %%%
data = []; %%%
if strcmp(monit_type, 'multi')   %%%
    W = zeros(1,n);
    static_set = zeros(1,n);
end
if strcmp(monit_type, 'max')
    W = zeros(p,n);
    static_set = zeros(p,n);
end
doNum0_set = cell(1,n);  %%%
doNum_set = cell(1,n);
for t = 1:w
    if strcmp(inter_type,'NI')
       j = 0;
       value = 0;
    end
    if strcmp(inter_type,'RI')|| strcmp(inter_type, 'AI')
        j = randi([0,p]);
        if j == 0
            value = 0;
        else
            value = doValue(j);
        end
    end
    if t < tau 
        data = [data; generate_data(A, mu, Sigma, 1, j, value)];
    else
        data = [data; generate_data(A_oc, mu_oc, Sigma_oc,1, j, value)];
    end
    doSet(t) = j;
end

mean_hat_set = cell(1,n);  %%%
cov_hat_set = cell(1,n);
mean_hat0_set = cell(1,n);
cov_hat0_set = cell(1,n);

for t = (w+1):n
    mean_hat = cell(p,1);   %%%
    cov_hat = cell(p,1);
    mean_hat0 = 0;
    cov_hat0 = 0;
    doNum = zeros(1,p);  %%%
    for l = 0:p 
        myseq = doSet((t-w):(t-1));
        index = find(myseq == l) + (t-w-1);
        temp_data = data(index,:);  
        if l == 0
            [mu_X,Sigma_X] = CalDisX(A, mu, Sigma, l, 0);
            temp_data = (temp_data - mu_X') * (eye(p) - A)' * inv(sqrtm(Sigma))'; %%
        else
            [mu_X,Sigma_X] = CalDisX(A, mu, Sigma, l, doValue(l));
            [A_do, mu_do,Sigma_do] = DoOperator(A,mu,Sigma,l,doValue(l));
            temp_data = (temp_data - mu_X') * (eye(p) - A_do)' * inv(sqrtm(Sigma))'; %%
        end
       


        if l == 0
            doNum0 = length(index);
        else
            doNum(l) = length(index);
        end

        if l == 0 && doNum0 >= 2        %%
            mean_hat0 = mean(temp_data);
            cov_hat0 = cov(temp_data,1);
        end
        if l ~=0 && doNum(l) >= 2 %%
            mean_hat{l} = mean(temp_data);
            cov_hat{l} = cov(temp_data,1);
        end
    end
    mean_hat_set{t} = mean_hat;
    cov_hat_set{t} = cov_hat;
    mean_hat0_set{t} = mean_hat0;
    cov_hat0_set{t} = cov_hat0;
    doNum_set{t} = doNum;
    doNum0_set{t} = doNum0;
    %%%%
    if strcmp(inter_type, 'NI')
        j = 0;
        value = 0;
    end
    if strcmp(inter_type, 'RI')
        j = randi([0,p]);
        if j == 0
            value = 0;
        else
            value = doValue(j);
        end
    end
    if strcmp(inter_type, 'AI')
        if seqN(t) == 1   %%%
            j = randi([0,p]);
            if j == 0
                value = 0;
            else
                value = doValue(j);
            end
        else   %%%%
            if strcmp(monit_type, 'multi')
                myKLdis0 = -1;
                myKLdis = -1*ones(1,p);
                if doNum0 >= p   %%%
                    myKLdis0 = normalKL(mean_hat0, cov_hat0);
                end
                for l = 1:p
                    if doNum(l) >= (p-1) %%%
                        temp_mean = mean_hat{l};
                        temp_mean(l) = [];
                        temp_cov = cov_hat{l};
                        temp_cov(l,:) = [];
                        temp_cov(:,l) = [];
                        myKLdis(l) = normalKL(temp_mean,temp_cov);
                    end
                end
                [Mvalue, Mindex] = max(myKLdis);
                if Mvalue >= myKLdis0
                    j = Mindex;
                    value = doValue(j);
                else
                    j = 0;
                    value = 0;
                end

            end
            if strcmp(monit_type, 'max')
                myKLdis0 = -1*ones(1,p);
                myKLdis  = -1 * ones(p,p); 
                if doNum0 >=2 
                    for l = 1:p
                        myKLdis0(l) = normalKL(mean_hat0(l), cov_hat0(l,l));
                    end
                end
                for l = 1:p
                    if doNum(l) >= 2
                        E = inv(eye(p) - A); %% 。
                        all_set = 1:p;
                        anc_set = find(E(l,1:(l-1))); 
                        no_anc_set = setdiff(all_set, [anc_set,l]); 
                        if isempty(no_anc_set)  %%%%
                            continue
                        end
                        for k = no_anc_set
                            myKLdis(l,k) = normalKL(mean_hat{l}(k), cov_hat{l}(k,k));
                        end
                    end
                end
                vec_myKLdis = max(myKLdis,[],2); %%
                vec_myKLdis0 = max(myKLdis0);
                [Mvalue,Mindex] = max(vec_myKLdis);
                if Mvalue >= vec_myKLdis0
                    j = Mindex;
                    value = doValue(j);
                else
                    j = 0;
                    value = 0;
                end
            end
        end
    end
    doSet(t) = j;
    %%%%
    if t < tau 
        data = [data; generate_data(A, mu, Sigma, 1, j, value)];
    else
        data = [data; generate_data(A_oc, mu_oc, Sigma_oc,1, j, value)];
    end
    if  strcmp(monit_type, 'multi')
        static = 0;
        if doSet(t) == 0 && doNum0 >= 2%p
            temp_cov = cov_hat0;
            [~,nonpositi] = chol(temp_cov);  %%
        end
        if doSet(t) == 0 &&doNum0 >= p && ~nonpositi  %%
            temp_mean = mean_hat0;
            %if doNum0 >= p && ~nonpositi
                temp_cov = cov_hat0;
            %else
            %    temp_cov = cov_hat0 + 0.01 *eye(p);
            %end
            temp_data = data(t,:);
            [mu_X,Sigma_X] = CalDisX(A, mu, Sigma, doSet(t), 0);
            temp_data = (temp_data - mu_X') * (eye(p) - A)' * inv(sqrtm(Sigma))'; 
            static = log(mvnpdf(temp_data, temp_mean,temp_cov)) - log(mvnpdf(temp_data, zeros(1,p),eye(p)));
        end
        
        if doSet(t) ~=0 && doNum(doSet(t)) >= 2 %(p-1)
             j = doSet(t);
             temp_cov = cov_hat{j};
             temp_cov(j,:) = [];
             temp_cov(:,j) = [];
             [~,nonpositi] = chol(temp_cov);  %%
        end
        if doSet(t) ~=0 && doNum(doSet(t)) >= (p-1)  && ~nonpositi %%
            j = doSet(t);
            temp_mean = mean_hat{j};
            temp_mean(j) = [];
            %if doNum(j) >= (p-1) && ~nonpositi
                temp_cov = cov_hat{j};
                temp_cov(j,:) = [];
                temp_cov(:,j) = [];
            %else
            %    temp_cov = cov_hat{j};
            %    temp_cov(j,:) = [];
            %    temp_cov(:,j) = [];
            %    temp_cov = temp_cov + 0.01*eye(p-1);
            %end
            
            temp_data = data(t,:);
            [mu_X,Sigma_X] = CalDisX(A, mu, Sigma, j, doValue(j));
            [A_do, mu_do,Sigma_do] = DoOperator(A,mu,Sigma,j,doValue(j));
            temp_data = (temp_data - mu_X') * (eye(p) - A_do)' * inv(sqrtm(Sigma))';
            temp_data(j) = [];
            static = log(mvnpdf(temp_data, temp_mean,temp_cov)) - log(mvnpdf(temp_data, zeros(1,p-1),eye(p-1)));
        end
        

        static_set(t) = static;
        W(t) = max([W(t-1),0]) + static;
        if W(t) >  b  
            T = t;  
            break
        end
    end
    if strcmp(monit_type, 'max')
        static = zeros(1,p);
        if doSet(t) == 0 && doNum0 >= 2  %%%
            temp_mean = mean_hat0;
            temp_cov = cov_hat0;
            temp_data = data(t,:);
            [mu_X,Sigma_X] = CalDisX(A, mu, Sigma, doSet(t), 0);
            temp_data = (temp_data - mu_X') * (eye(p) - A)' * inv(sqrtm(Sigma))'; 
        end
        if doSet(t) ~=0 && doNum(doSet(t)) >=2 %%
            j = doSet(t);
            temp_mean = mean_hat{j};
            temp_cov = cov_hat{j};
            temp_data = data(t,:);
            [mu_X,Sigma_X] = CalDisX(A, mu, Sigma, j, doValue(j));
            [A_do, mu_do,Sigma_do] = DoOperator(A,mu,Sigma,j,doValue(j));
            temp_data = (temp_data - mu_X') * (eye(p) - A_do)' * inv(sqrtm(Sigma))';
        end
        if (doSet(t) == 0 && doNum0 >= 2) || (doSet(t) ~=0 && doNum(doSet(t)) >=2)
            for l = 1:p
                if l == doSet(t)
                   continue 
                end
                static(l) = log(normpdf(temp_data(l),temp_mean(l),temp_cov(l,l)))- log(normpdf(temp_data(l),0,1));
            end
        end
        static_set(:,t) = static';
        for l = 1:p
            W(l,t) = max([W(l,t-1),0]) + static(l);
        end
        if max(W(:,t)) > b %%%
            T = t;
            break
        end
    end

end

if strcmp(monit_type, 'max')
    W0 = W;
    W = max(W0,[],1);
end

 fid = fopen('tt.txt','a'); %% count txt.
 fprintf(fid, '%g\n', tt);
 fclose(fid);
 return;
 end







