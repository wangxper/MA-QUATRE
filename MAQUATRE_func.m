function [gbest,gbestval,RecordT] = MAQUATRE_func(fhd,D,ps,iter_max,Xmin,Xmax,varargin)
%% MAQUATRE_FUNC 
    % random seeds
    stm = RandStream('swb2712','Seed',sum(100*clock));
    RandStream.setGlobalStream(stm);
    % arguments passing
    gen = iter_max; 
    % ps*D matrix
    if length(Xmin) ==1
        Rmin = repmat(Xmin,1,D);
        Rmax = repmat(Xmax,1,D);
    end
    % Boundary
    VRmin = repmat(Rmin,ps,1);
    VRmax = repmat(Rmax,ps,1);
    
    fidvec = cell2mat(varargin);
    fid = fidvec(1);
    runid = fidvec(2);   
    %global best
    targetbest = [300;400;600;800;900;1800;2000;2200;2300;2400;2600;2700];

    F = 0.7;
    R_div = 0.5; % Out Rate
    RP = 0.5;    % Random Proportion in out-communication
    name = ['MAQUATRE_fid_',num2str(fid),'_',num2str(D),'D','.dat'];
    fout = fopen(name,'a');  
%% set plot-points
  if  runid ==1
        for i=0:50 % total 51 points
            if i==0
                iteration=0;
                fprintf(fout,'%s:%s\t','iteration_F',num2str(fid));
            else
                iteration=gen/50*i;
            end
            fprintf(fout,'%d\t',iteration);
        end       
        fprintf(fout,'\n');
  end
%%    
    tic;
    % initialize   
    pos = VRmin + (VRmax-VRmin).*rand(ps,D);
    % pos = randn(ps,D);
    eval = feval(fhd,pos',varargin{:}); % evaluation

    index = ceil(ps/2);
    p1pos  = pos(1:index,:);  % elementary membrane 1
    p1eval = eval(1:index);
    p2pos  = pos(index+1:ps,:);  % elementary membrane 2
    p2eval = eval(index+1:ps); 
    
    [gbestval,gbestid] = min(eval);
    gbest = pos(gbestid,:);
    
      fprintf(fout,'%s\t%.15f\t',num2str(runid),gbestval-targetbest(fid));

    for iter = 2:gen
       %% Out-communication rule:objects divide from membranes 1 and 2
        % membrane 1 out-operation
        [p1eval,index1] = sort(p1eval); % sort
        p1pos= p1pos(index1,:); 
        
        N_p1 = length(p1eval);  % objects number in membrane 1 
        N_p1_d = ceil(R_div * N_p1);   % objects out number in membrane 1
        
        ri1=randperm(N_p1,ceil(N_p1_d*RP));  % The proportion rate of output operation in the membrane 1, randomly selected       
        p1pos_r = p1pos(ri1,:);   % Randomly output objects in membrane 1
        p1eval_r = p1eval(1,ri1); % Fitness values of randomly output objects in membrane 1
        
        p1_index= 1:length(p1eval);
        p1_index(ismember(p1_index, ri1)) = [];
        
        wi1 = p1_index(length(p1_index)-(N_p1_d-length(ri1))+1:length(p1_index)); % Index of the output poorer solutions in membrane 1
        p1pos_w = p1pos(wi1,:);   % the output poorer objects in membrane 1
        p1eval_w = p1eval(1,wi1); % Fitness values of the output poorer objects in membrane 1
        
        p1_index(ismember(p1_index, wi1)) = [];       
        p1pos= p1pos(p1_index,:);   % Objects contained after membrane 1 out-communication rule
        p1eval= p1eval(1,p1_index); % Fitness values of objects contained after membrane 1 out-communication rule
        
       % membrane 2 out-operation - Refer to the annotation of membrane 1
        [p2eval,index2] = sort(p2eval); % sort
        p2pos= p2pos(index2,:); 
        
        N_p2 = length(p2eval);  
        N_p2_d = ceil(R_div * N_p2);   
        
        ri2=randperm(N_p2,ceil(N_p2_d*RP));      
        p2pos_r = p2pos(ri2,:);   
        p2eval_r = p2eval(1,ri2); 
        
        p2_index= 1:length(p2eval);
        p2_index(ismember(p2_index, ri2)) = [];
        
        wi2 = p2_index(length(p2_index)-(N_p2_d-length(ri2))+1:length(p2_index)); 
        p2pos_w = p2pos(wi2,:);   
        p2eval_w = p2eval(1,wi2); 
        
        p2_index(ismember(p2_index, wi2)) = [];       
        p2pos= p2pos(p2_index,:);   
        p2eval= p2eval(1,p2_index); 
        % membrane 3  skin membrane
        p3pos= [p1pos_r;p1pos_w;p2pos_r;p2pos_w];
        p3eval= [p1eval_r,p1eval_w,p2eval_r,p2eval_w];       

       %% Generating evolution matrix M.        
        K1 = floor(length(p1eval)/D);
        K2 = floor(length(p2eval)/D);
        K3 = floor(length(p3eval)/D);
        
        % p1_M membrane 1   M matrix
        p1_M = zeros(length(p1eval),D);
        for j = 0:K1-1
            p1_M(j*D+1:j*D+D,:)=tril(ones(D,D),0);
        end
        % ps1/D Assignment when not divided by integer
           p1_M(K1*D+1:length(p1eval),:) = tril(ones(length(p1eval)-K1*D,D),0);
        
        % permutation of the matrix
        for t = 1:length(p1eval)
            p1_M(t,:) = p1_M(t,randperm(D));
        end
        p1_M = p1_M(randperm(size(p1_M,1))',:);
        
        % p2_M membrane 2   M matrix
        p2_M = zeros(length(p2eval),D);
        for j = 0:K2-1
            p2_M(j*D+1:j*D+D,:)=tril(ones(D,D),0);
        end
        % ps2/D Assignment when not divided by integer
           p2_M(K2*D+1:length(p2eval),:) = tril(ones(length(p2eval)-K2*D,D),0);
        
        % permutation of the matrix
        for t = 1:length(p2eval)
            p2_M(t,:) = p2_M(t,randperm(D));
        end
        p2_M = p2_M(randperm(size(p2_M,1))',:);
        
        % p3_M  skin membrane 3   M matrix
        p3_M = zeros(length(p3eval),D);
        for j = 0:K3-1
            p3_M(j*D+1:j*D+D,:)=tril(ones(D,D),0);
        end
        % ps3/D Assignment when not divided by integer
           p3_M(K3*D+1:length(p3eval),:) = tril(ones(length(p3eval)-K3*D,D),0);
        
        % permutation of the matrix
        for t = 1:length(p3eval)
            p3_M(t,:) = p3_M(t,randperm(D));
        end
        p3_M = p3_M(randperm(size(p3_M,1))',:);
       %% Evolution Operator 
        % p1Rmin(ps1*D) Boundary of membrane 1
        p1Rmin = repmat(Rmin,length(p1eval),1);
        p1Rmax = repmat(Rmax,length(p1eval),1);
        
        p1pos_r1 = p1pos(randperm(length(p1eval))',:);
        p1pos_r2 = p1pos(randperm(length(p1eval))',:);
        p1pos_r3 = p1pos(randperm(length(p1eval))',:);
        
        p1pos_new = p1pos_r1 + F*(p1pos_r2-p1pos_r3);   % rand/1
        
        p1pos_new(~p1_M) = p1pos(~p1_M);
        % Boundary Limits
        p1pos_new = ((p1pos_new>=p1Rmin)&(p1pos_new<=p1Rmax)).*p1pos_new+(p1pos_new<p1Rmin).*(p1Rmin+p1pos)/2+(p1pos_new>p1Rmax).*(p1Rmax+p1pos)/2;
         
        p1eval_new = feval(fhd,p1pos_new',varargin{:});
        
        bin1 = (p1eval > p1eval_new)'; 
        p1pos(bin1==1,:) = p1pos_new(bin1==1,:);
        p1eval(bin1==1) = p1eval_new(bin1==1);
        
        [p1_gbestval,p1_gbestid] = min(p1eval);
        p1_gbest = p1pos(p1_gbestid,:);
        if p1_gbestval < gbestval  
            gbestval = p1_gbestval;
            gbest = p1_gbest;      
        end
       
        % p2Rmin(ps2*D) Boundary of membrane 2
        p2Rmin = repmat(Rmin,length(p2eval),1);
        p2Rmax = repmat(Rmax,length(p2eval),1);
        
        [~,p2_gbestid] = min(p2eval);
        p2gbest = p2pos(p2_gbestid,:);
        p2gbestrep= repmat(p2gbest,length(p2eval),1);   
        
        p2pos_r1 = p2pos(randperm(length(p2eval))',:);
        p2pos_r2 = p2pos(randperm(length(p2eval))',:);
        
        p2pos_new = p2gbestrep + F*(p2pos_r1-p2pos_r2); % best/1
              
        p2pos_new(~p2_M) = p2pos(~p2_M);
        
        p2pos_new = ((p2pos_new>=p2Rmin)&(p2pos_new<=p2Rmax)).*p2pos_new+(p2pos_new<p2Rmin).*(p2Rmin+p2pos)/2+(p2pos_new>p2Rmax).*(p2Rmax+p2pos)/2;
       
        p2eval_new = feval(fhd,p2pos_new',varargin{:});
        
        bin2 = (p2eval > p2eval_new)'; 
        p2pos(bin2==1,:) = p2pos_new(bin2==1,:);
        p2eval(bin2==1) = p2eval_new(bin2==1);
        
        [p2_gbestval,p2_gbestid] = min(p2eval);
        p2_gbest = p2pos(p2_gbestid,:);   
        if p2_gbestval < gbestval  
            gbestval = p2_gbestval;
            gbest = p2_gbest;      
        end
        
        % p3Rmin(ps3*D) Boundary of skin membrane 3
        p3Rmin = repmat(Rmin,length(p3eval),1);
        p3Rmax = repmat(Rmax,length(p3eval),1);
        
        [~,p3_gbestid] = min(p3eval);
        p3gbest = p3pos(p3_gbestid,:);
        p3gbestrep= repmat(p3gbest,length(p3eval),1);   
        
        p3pos_r1 = p3pos(randperm(length(p3eval))',:);
        p3pos_r2 = p3pos(randperm(length(p3eval))',:);
        
        p3pos_new = p3pos + F*(p3gbestrep-p3pos) + F*(p3pos_r1-p3pos_r2); % target to best/1
        
        p3pos_new(~p3_M) = p3pos(~p3_M);
        
        p3pos_new = ((p3pos_new>=p3Rmin)&(p3pos_new<=p3Rmax)).*p3pos_new+(p3pos_new<p3Rmin).*(p3Rmin+p3pos)/2+(p3pos_new>p3Rmax).*(p3Rmax+p3pos)/2;
        
        p3eval_new = feval(fhd,p3pos_new',varargin{:});
        
        bin3 = (p3eval > p3eval_new)'; 
        p3pos(bin3==1,:) = p3pos_new(bin3==1,:);
        p3eval(bin3==1) = p3eval_new(bin3==1);
        
        [p3_gbestval,p3_gbestid] = min(p3eval);
        if p3_gbestval < gbestval  
            gbestval = p3_gbestval;
            gbest = p3pos(p3_gbestid,:);      
        end

       %% In-communication rule
       % Calculate the distance from the updated out object to the center of each membrane to membranes 1 and 2
        DTable = CDistance(p1pos,p2pos,p3pos);    
        p1_indices = find(DTable == 1);
        p2_indices = find(DTable == 2);
         
        p1pos_f = p3pos(p1_indices,:);   % membrane 1: In-operation objects from the skin membrane 3  
        p1eval_f = p3eval(1,p1_indices); % membrane 1: Objects fitness of In-operation objects
        
        p2pos_f = p3pos(p2_indices,:);   % membrane 2: In-operation objects from the skin membrane 3  
        p2eval_f = p3eval(1,p2_indices); % membrane 2: Objects fitness of In-operation objects
              
        p1pos= [p1pos;p1pos_f];        % Object of membrane 1 after In-communication rule
        p1eval= [p1eval,p1eval_f];    
        
        p2pos= [p2pos;p2pos_f];        % Object of membrane 2 after In-communication rule
        p2eval= [p2eval,p2eval_f];   
      %% After In-communication rule, there should be at least 4 solutions in membranes 1 and 2. 
       % If there are less than 4, the other membrane should be divided in half randomly; 
       % setting 4 solutions is to ensure that membranes 1 and 2 can be divided again.    
        if length(p1eval)<4
            
           p2_add=randperm(length(p2eval),ceil(length(p2eval)/2));  % Randomly selected objects from membrane 2 are added to membrane 1   
           p1pos_add = p2pos(p2_add,:);   % Additional objects
           p1eval_add = p2eval(1,p2_add); % Additional objects fitness 
           
           p1pos= [p1pos;p1pos_add];      % Objects of membrane 1 after addition
           p1eval= [p1eval,p1eval_add];   
          
           p2pos(p2_add,:)= [];           % Objects of membrane 2 after addition
           p2eval(p2_add)= [];  
           
        elseif length(p2eval)<4
           p1_add=randperm(length(p1eval),ceil(length(p1eval)/2));  % Randomly selected objects from membrane 1 are added to membrane 2  
           p2pos_add = p1pos(p1_add,:);   % Additional objects
           p2eval_add = p1eval(1,p1_add); % Additional objects fitness 

           p1pos(p1_add,:)= [];           % Objects of membrane 1 after addition
           p1eval(p1_add)= [];  
           
           p2pos= [p2pos;p2pos_add];      % Objects of membrane 2 after addition
           p2eval= [p2eval,p2eval_add];   
            
        end
           
        if mod(iter,gen/50)==0
             fprintf(fout,'%.15f\t',gbestval-targetbest(fid));  
        end
        
    end
    RecordT = toc;
    
    fprintf(fout,'\n');
    fclose(fout);
end

function [DTable]=CDistance(p1pos,p2pos,p3pos)
DTable=zeros(size(p3pos,1),1);
Tabled=zeros(size(p3pos,1),2);
p1_center=mean(p1pos);
p2_center=mean(p2pos);

 for i=1:size(p3pos,1)
     for j=1:2
         if j==1
         Tabled(i,j)=norm(p3pos(i,:)-p1_center); % table distance
         else
         Tabled(i,j)=norm(p3pos(i,:)-p2_center); % table distance
         end
     end
        [~,index] = sort(Tabled(i,:));
        DTable(i,1)=index(1,1);    
 end
end 