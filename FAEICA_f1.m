tic;
% FAEICA for f1 (Sphere)
clear;
close all;
%MaxRun = 30;
MaxRun = 1;
%MaxFEs = 20000;
MaxFEs = 840;
FEsRuns = zeros(MaxRun, 1);
CostRuns=zeros(MaxRun, 1);

for nRun=1:MaxRun
disp(['Independant Run : ' num2str(nRun)]);
%% Problem Definition
CostFunction = @f1; % Sphere
nVar = 10;
down=[-100]; %#ok<*NBRAK>
up=[+100];
VarSize = [1 nVar];

% Initial settings for FAEICA
MaxIt = 9472; % Stopping criterion: FEs at each iteration*MaxIt == 20,000 for f1, so MaxIt can be any number larger than 20,000/FEs at each iteration
Alpha = 10; % This Alpha is used for the change range in the velocity component of each solution
nPop = 100;
nEmp = 20;
nCol = nPop - nEmp;

C=0.3;
% Probability of applying Velocity Divergence
PVD = 1;

% Probability of applying Velocity Adaptation
PVA = 1;

% Probability of applying Information Sharing between empires
PLS = 1;

for j = 1:nVar
  GlobalBest.Position(j) = rand.*(up-down) + down;
end
GlobalBest.Velocity = zeros(VarSize);
GlobalBest.Cost = Inf;
GlobalBest.Power = 0; % The worst value

%% Initialization
empty_country.Position=[];
empty_country.Velocity=[];
empty_country.Cost=[];
empty_country.Power=[];

empty_country.Best.Position=[];
empty_country.Best.Velocity=[];
empty_country.Best.Cost=[];
empty_country.Best.Power=[];

country=repmat(empty_country,nPop,1);

BestPower = zeros(MaxIt,1);
BestCost = zeros(MaxIt,1);
GlobalBestCost = zeros(1, 1);
GlobalBestPosition = zeros(1, nVar);
convergence1 = zeros(MaxIt, 1);
    
  FEs = 0;
  for i=1:nPop
      
      for j = 1:nVar
          country(i).Position(j) = rand.*(up-down) + down;
      end
      
      country(i).Velocity = zeros(VarSize);
      [country(i).Cost, country(i).Power] = CostFunction(country(i).Position);      
      country(i).Best.Position = country(i).Position;
      country(i).Best.Velocity = country(i).Velocity;                 
      country(i).Best.Cost = country(i).Cost;      
      country(i).Best.Power = country(i).Power;            
      
      if country(i).Best.Power > GlobalBest.Power
              GlobalBest = country(i).Best;
      end                                           
  end
    
%%  Form imperialists and colonies    
    % Sort countries
    [~,index]=sort([country.Power],'descend');
    country=country(index);
    % Assign of Colonies and Imperialists
    imp=country(1:nEmp);
    col=country(nEmp+1:end);

    empty_empire.Imp=[];
    empty_empire.Col=repmat(empty_country,0,1);
    empty_empire.nCol=0;
    emp=repmat(empty_empire,nEmp,1);

   % Assign Imperialists
        for k=1:nEmp
            emp(k).Imp=imp(k);
        end

    % Assign Colonies
        Np = [imp.Power]./max([imp.Power]);
        SNp = sum(Np);
        Npp = Np./SNp;

    for Num=1:nCol
         k = RouletteWheelSelection(Npp);
         check = isempty(emp(k).Col);
         if check == 1
             emp(k).Col = col(Num);
         else
            emp(k).Col = [emp(k).Col; col(Num)];
         end
         emp(k).nCol = emp(k).nCol + 1;                 
    end
        
    index = find([emp.nCol]==0);
    if index>=1
        for tedademp=1:numel(index)
            [quan, ind]=max([emp.nCol]);
            randcol = randi(quan);
            emp(index(tedademp)).Col = emp(ind).Col(randcol);
            emp(index(tedademp)).nCol = 1;
            emp(ind).Col(randcol) = [];
            emp(ind).nCol = emp(ind).nCol - 1;             
        end  
    end
    
 %% Main Loop   
    for it=1:MaxIt
    fprintf('Iteration time is %d \n', it)
% FIS1 is used by Global Learning, Universal global best diversity,
% Differential evolutionary-based local search to adapt the parameters at
% each iteration time

% FIS2 is used for dynamic selection of operators at each time window
if mod(it, MaxFEs/10) == 0 % tw = MaxFEs/10;
%% Fuzzy Adaptive Operator Selection
    BestPower_tw = BestPower(it-10+1:it-1);
    Delta_BestPower_tw =  max(BestPower_tw) - min(BestPower_tw);
    Stagnation = 1 - (Delta_BestPower_tw / max(BestPower_tw));
    Stagnation = min(Stagnation, 1);
    UFuzzy = [Stagnation; PVA; PVD; PLS];
    FISMAT = readfis('FAOS.fis');
    Y = evalfis(FISMAT, UFuzzy);
    PVA = Y(1,1);
    PVD = Y(1,2);
    PLS = Y(1,3);
end

%% Global learning?based velocity adaptation
    if rand <= PVA
    for k = 1:numel(emp)
        for col = 1:numel(emp(k).Col)
            NP = (abs(emp(k).Imp.Power - emp(k).Col(col).Power))/emp(k).Imp.Power;
            NP = min(NP,1);
            NFEs = FEs/MaxFEs;            
            UFuzzy = [NP; NFEs];
            FISMAT = readfis('FAGLVA.fis');
            Y = evalfis(FISMAT, UFuzzy);
            % Social learning parameters {w, c2, Beta}
            % Cognitive learning parameters {c1}
            w = Y(1,1);
            c1 = Y(1,2);
            Beta = Y(1,3);
            c2 = Y(1,4);
            
            tempCol1 = emp(k).Col(col);
            tempImp1 = emp(k).Imp;                                   
            emp(k).Col(col).Velocity = (w.*emp(k).Col(col).Velocity)+(Beta.*rand(VarSize)).*(emp(k).Imp.Position-emp(k).Col(col).Position) + (c1.*rand(VarSize)).*(emp(k).Col(col).Best.Position-emp(k).Col(col).Position) + (c2.*rand(VarSize)).*(GlobalBest.Position-emp(k).Col(col).Position);                                              
            
            [VelMin, VelMax] = VelLimit(GlobalBest.Position, emp(k).Col(col).Position, it, up, down, Alpha);
            emp(k).Col(col).Velocity = min(max(emp(k).Col(col).Velocity,VelMin),VelMax);
            [VelMin, VelMax] = VelLimit(GlobalBest.Position, emp(k).Imp.Position, it, up, down, Alpha);
            emp(k).Imp.Velocity = min(max(emp(k).Imp.Velocity,VelMin),VelMax);                   
            emp(k).Col(col).Position = emp(k).Col(col).Position + emp(k).Col(col).Velocity;
            emp(k).Imp.Position = emp(k).Imp.Position + emp(k).Imp.Velocity;
            for flg=1:nVar
                if emp(k).Col(col).Position(flg) < down || emp(k).Col(col).Position(flg) > up 
                    emp(k).Col(col).Velocity(flg) = -emp(k).Col(col).Velocity(flg);              
                end
            end
            for flg=1:nVar
                if emp(k).Imp.Position(flg) < down || emp(k).Imp.Position(flg) > up %#ok<*BDSCI>
                    emp(k).Imp.Velocity(flg) = -emp(k).Imp.Velocity(flg);              
                end
            end
            for flg=1:nVar
                    emp(k).Imp.Position(flg) = min(up,max(down,emp(k).Imp.Position(flg))); % Bound the new location
            end
            for flg=1:nVar
                    emp(k).Col(col).Position(flg) = min(up,max(down,emp(k).Col(col).Position(flg))); % Bound the new location
            end
          
            [emp(k).Col(col).Cost, emp(k).Col(col).Power] = CostFunction(emp(k).Col(col).Position);
            FEs = FEs + 1;            
            
            if emp(k).Col(col).Power >= emp(k).Col(col).Best.Power
                emp(k).Col(col).Best.Position = emp(k).Col(col).Position;
                emp(k).Col(col).Best.Velocity = emp(k).Col(col).Velocity;
                emp(k).Col(col).Best.Cost = emp(k).Col(col).Cost;
                emp(k).Col(col).Best.Power = emp(k).Col(col).Power;
            end
            
            if emp(k).Col(col).Power >= GlobalBest.Power
                temp = GlobalBest;
                GlobalBest.Position = emp(k).Col(col).Position;
                GlobalBest.Velocity = emp(k).Col(col).Velocity;
                GlobalBest.Cost = emp(k).Col(col).Cost;
                GlobalBest.Power = emp(k).Col(col).Power;
                emp(k).Col(col).Position = temp.Position;
                emp(k).Col(col).Velocity = temp.Velocity;
                emp(k).Col(col).Cost = temp.Cost;
                emp(k).Col(col).Power = temp.Power;
            end

            % Exchange a colony with its local imperialist if the colony has more power than its local imperialist
            if emp(k).Col(col).Power > emp(k).Imp.Power
                [emp(k).Imp, emp(k).Col(col)] = deal(emp(k).Col(col), emp(k).Imp);                
            end
            
            [emp(k).Imp.Cost, emp(k).Imp.Power] = CostFunction(emp(k).Imp.Position);
            FEs = FEs + 1;
            

            if emp(k).Imp.Power >= emp(k).Imp.Best.Power
                emp(k).Imp.Best.Position = emp(k).Imp.Position;
                emp(k).Imp.Best.Velocity = emp(k).Imp.Velocity;
                emp(k).Imp.Best.Cost = emp(k).Imp.Cost;
                emp(k).Imp.Best.Power = emp(k).Imp.Power;
            end    
            
            if emp(k).Imp.Power >= GlobalBest.Power
                temp = GlobalBest;
                GlobalBest.Position = emp(k).Imp.Position;
                GlobalBest.Velocity = emp(k).Imp.Velocity;
                GlobalBest.Cost = emp(k).Imp.Cost;
                GlobalBest.Power = emp(k).Imp.Power;
                emp(k).Imp.Position = temp.Position;
                emp(k).Imp.Velocity = temp.Velocity;
                emp(k).Imp.Cost = temp.Cost;
                emp(k).Imp.Power = temp.Power;
            end        
        end
    end
    end
            
%% Universal global best diversity
    if rand <= PVD
    for k = 1:numel(emp)
        for col = 1:numel(emp(k).Col)
            NP = abs(emp(k).Imp.Power - emp(k).Col(col).Power)/emp(k).Imp.Power;
            NP = min(NP,1);
            NFEs = FEs/MaxFEs;             
            UFuzzy = [NP; NFEs];
            FISMAT = readfis('FAUDVD.fis');
            Y = evalfis(FISMAT, UFuzzy);
            Pdiv = Y(1,1);
            Np = [imp.Power]./max([imp.Power]);
            SNp = sum(Np);     
            if SNp == 0
                Npp = Np;
            else
            Npp = Np./SNp;
            end             

            tempCol = emp(k).Col(col);
            for d = 1:nVar
                randk = RouletteWheelSelection(Npp);
                while numel(emp(randk).Col) == 0
                    randk = RouletteWheelSelection(Npp);
                end
                Npcol = [emp(randk).Col.Power]./max([emp(randk).Col.Power]);
                SNpcol = sum(Npcol);
                if SNpcol == 0
                    Nppcol = Npcol;
                else
                Nppcol = Npcol./SNpcol;            
                end
                status = isnan(Nppcol);
                if (Nppcol ~= 0) & (status ~= 1) %#ok<*AND2>
                    randcol = RouletteWheelSelection(Nppcol);
                else
                    randcol = randi(numel(emp(randk).Col));
                end
                emp(k).Col(col).Velocity(d) = emp(k).Col(col).Velocity(d) + rand.*(emp(randk).Col(randcol).Best.Position(d) - emp(k).Col(col).Position(d));
            end

            tempGL = GlobalBest;
            for d = 1:nVar
                randk = RouletteWheelSelection(Npp);
                if rand < PVD
                    GlobalBest.Velocity(d) = GlobalBest.Velocity(d) + rand.*(emp(randk).Imp.Best.Position(d) - GlobalBest.Position(d));
                end
            end
                                  
            % Constraints are applied to the velocities
            [VelMin, VelMax] = VelLimit(GlobalBest.Position, emp(k).Col(col).Position, it, up, down, Alpha);
            emp(k).Col(col).Velocity = min(max(emp(k).Col(col).Velocity,VelMin),VelMax);           
            [VelMin, VelMax] = VelLimit(GlobalBest.Position, emp(k).Imp.Position, it, up, down, Alpha);
            emp(k).Imp.Velocity = min(max(emp(k).Imp.Velocity,VelMin),VelMax);
            VelMax = 0; VelMin = 0;           
            
            GlobalBest.Velocity = min(max(GlobalBest.Velocity,VelMin),VelMax);

            emp(k).Col(col).Position = emp(k).Col(col).Position + emp(k).Col(col).Velocity;
            emp(k).Imp.Position = emp(k).Imp.Position + emp(k).Imp.Velocity;
            GlobalBest.Position = GlobalBest.Position + GlobalBest.Velocity;
            
            for flg=1:nVar
                if emp(k).Col(col).Position(flg) < down | emp(k).Col(col).Position(flg) > up %#ok<*OR2>
                    emp(k).Col(col).Velocity(flg) = -emp(k).Col(col).Velocity(flg);              
                end
            end
            for flg=1:nVar
                    emp(k).Col(col).Position(flg) = min(up,max(down,emp(k).Col(col).Position(flg))); % Bound the new location
            end
            
            for flg=1:nVar
                if emp(k).Imp.Position(flg) < down | emp(k).Imp.Position(flg) > up
                    emp(k).Imp.Velocity(flg) = -emp(k).Imp.Velocity(flg);              
                end
            end
            for flg=1:nVar
                    emp(k).Imp.Position(flg) = min(up,max(down,emp(k).Imp.Position(flg))); % Bound the new location
            end
            
            for flg=1:nVar
                if GlobalBest.Position(flg) < down | GlobalBest.Position(flg) > up
                    GlobalBest.Velocity(flg) = -GlobalBest.Velocity(flg);              
                end
            end
            for flg=1:nVar
                    GlobalBest.Position(flg) = min(up,max(down,GlobalBest.Position(flg))); % Bound the new location
            end
            [emp(k).Col(col).Cost, emp(k).Col(col).Power] = CostFunction(emp(k).Col(col).Position);
            FEs = FEs + 1;
                
            % Update personal best of colonies
            if emp(k).Col(col).Power >= emp(k).Col(col).Best.Power
                emp(k).Col(col).Best.Position = emp(k).Col(col).Position;
                emp(k).Col(col).Best.Velocity = emp(k).Col(col).Velocity;
                emp(k).Col(col).Best.Cost = emp(k).Col(col).Cost;
                emp(k).Col(col).Best.Power = emp(k).Col(col).Power;
            end

            % Update the Global Best imperialist using colonies
            if emp(k).Col(col).Power >= GlobalBest.Power
                temp = GlobalBest;
                GlobalBest.Position = emp(k).Col(col).Position;
                GlobalBest.Velocity = emp(k).Col(col).Velocity;
                GlobalBest.Cost = emp(k).Col(col).Cost;
                GlobalBest.Power = emp(k).Col(col).Power;
                emp(k).Col(col).Position = temp.Position;
                emp(k).Col(col).Velocity = temp.Velocity;
                emp(k).Col(col).Cost = temp.Cost;
                emp(k).Col(col).Power = temp.Power;
            end

            % Exchange a colony with its local imperialist if the colony has more power than its local imperialist
            if emp(k).Col(col).Power > emp(k).Imp.Power
                [emp(k).Imp, emp(k).Col(col)] = deal(emp(k).Col(col), emp(k).Imp);
            end        

            [emp(k).Imp.Cost, emp(k).Imp.Power] = CostFunction(emp(k).Imp.Position);
            FEs = FEs + 1;
            
           % Update personal best of the imperialist1
            if emp(k).Imp.Power >= emp(k).Imp.Best.Power
                emp(k).Imp.Best.Position = emp(k).Imp.Position;
                emp(k).Imp.Best.Velocity = emp(k).Imp.Velocity;
                emp(k).Imp.Best.Cost = emp(k).Imp.Cost;
                emp(k).Imp.Best.Power = emp(k).Imp.Power;
            end

            % Update the Global Best imperialist using colonies
            if emp(k).Imp.Power >= GlobalBest.Power
                temp = GlobalBest;
                GlobalBest.Position = emp(k).Imp.Position;
                GlobalBest.Velocity = emp(k).Imp.Velocity;
                GlobalBest.Cost = emp(k).Imp.Cost;
                GlobalBest.Power = emp(k).Imp.Power;
                emp(k).Imp.Position = temp.Position;
                emp(k).Imp.Velocity = temp.Velocity;
                emp(k).Imp.Cost = temp.Cost;
                emp(k).Imp.Power = temp.Power;
            end        
 
            [GlobalBest.Cost, GlobalBest.Power] = CostFunction(GlobalBest.Position);
            FEs = FEs + 1;            
            % Global Best is sensitive, selection phase is required!       
            if tempGL.Power > GlobalBest.Power
                 [GlobalBest, tempGL] = deal(tempGL, GlobalBest);                
            end                         
        end
    end   
    end
    
%% Differential evolutionary-based local search
   if rand <= PLS
        for k = 1:numel(emp)
            if numel(emp(k).Col) ~= 0
                    temp = [emp(k).Col.Power];
                    [~, indexwCol] = min(temp);
                    wCol = emp(k).Col(indexwCol);   
                    line_found = find([emp.nCol]~=0);
                    rand1 = line_found(randi(numel(find([emp.nCol]~=0))));

                    Colr1 = emp(rand1).Col(randi(numel(emp(rand1).Col)));
                        if Colr1.Position == wCol.Position
                            line_found = find([emp.nCol]~=0);
                            rand1 = line_found(randi(numel(find([emp.nCol]~=0))));                            
                            if numel(emp(rand1).Col) == 0
                                line_found = find([emp.nCol]~=0);
                                rand1 = line_found(randi(numel(find([emp.nCol]~=0))));
                            end
                        end
                    Colr1 = emp(rand1).Col(randi(numel(emp(rand1).Col)));
                        
                    line_found = find([emp.nCol]~=0);
                    rand2 = line_found(randi(numel(find([emp.nCol]~=0))));                   
                    Colr2 = emp(rand2).Col(randi(numel(emp(rand2).Col)));     
                    if Colr2.Position == wCol.Position
                        line_found = find([emp.nCol]~=0);
                        rand2 = line_found(randi(numel(find([emp.nCol]~=0))));
                        
                        if numel(emp(rand2).Col) == 0
                            line_found = find([emp.nCol]~=0);
                            rand2 = line_found(randi(numel(find([emp.nCol]~=0))));                            
                        end                            
                        Colr2 = emp(rand2).Col(randi(numel(emp(rand2).Col)));                        
                    end
                    
                    % The NP is calculated for the Colr3 and Imp
                    NP = abs(emp(k).Imp.Power - wCol.Power)/emp(k).Imp.Power;
                    NP = min(NP,1);
                    NFEs = FEs/MaxFEs;
                    UFuzzy = [NP; NFEs];
                    FISMAT = readfis('FADELS.fis');
                    Y = evalfis(FISMAT, UFuzzy);
                    F1 = Y(1,1);
                    F2 = Y(1,2);
                    pCR = Y(1,3);
                    MutantCol.Velocity = (rand*F1).*(emp(k).Imp.Position - wCol.Position) + (rand*F2).*(Colr1.Position - Colr2.Position);

                    % Constraints are applied to the velocities            
                    VelMax = +Alpha.*((up-down)./up); VelMin = -VelMax;                     
                    MutantCol.Velocity = min(max(MutantCol.Velocity,VelMin),VelMax);
                    
                    for d=1:nVar
                        if rand <= PLS
                            TrialCol.Velocity(d) = MutantCol.Velocity(d);
                        else
                            TrialCol.Velocity(d) = wCol.Velocity(d);
                        end
                    end
 
                    % Constraints are applied to the velocities            
                    VelMax = +Alpha.*((up-down)./up); VelMin = -VelMax;                     
                    TrialCol.Velocity = min(max(TrialCol.Velocity,VelMin),VelMax);

                    TrialCol.Position = wCol.Position + TrialCol.Velocity;
                    
                    for flg=1:nVar
                        if TrialCol.Position(flg) < down | TrialCol.Position(flg) > up
                            TrialCol.Velocity(flg) = -TrialCol.Velocity(flg);              
                        end
                    end
                    for flg=1:nVar
                            TrialCol.Position(flg) = min(up,max(down,TrialCol.Position(flg))); % Bound the new location
                    end                    
                    
                    [TrialCol.Cost, TrialCol.Power] = CostFunction(TrialCol.Position);
                    FEs = FEs + 1;
                    % Selection
                    if (TrialCol.Cost*(norm(emp(k).Imp.Position-TrialCol.Position))) < (wCol.Cost*(norm(emp(k).Imp.Position-wCol.Position)))
                        emp(k).Col(indexwCol).Position = TrialCol.Position;                    
                        emp(k).Col(indexwCol).Velocity = TrialCol.Velocity;
                        emp(k).Col(indexwCol).Cost = TrialCol.Cost;
                        emp(k).Col(indexwCol).Power = TrialCol.Power;
                    end

                    % Update personal best of the selected colony
                    if emp(k).Col(indexwCol).Power >= emp(k).Col(indexwCol).Best.Power
                        emp(k).Col(indexwCol).Best.Position = emp(k).Col(indexwCol).Position;
                        emp(k).Col(indexwCol).Best.Velocity = emp(k).Col(indexwCol).Velocity;
                        emp(k).Col(indexwCol).Best.Cost = emp(k).Col(indexwCol).Cost;
                        emp(k).Col(indexwCol).Best.Power = emp(k).Col(indexwCol).Power;
                    end

                    % Update the Global Best imperialist using colonies
                    if emp(k).Col(indexwCol).Power >= GlobalBest.Power
                        temp = GlobalBest;
                        GlobalBest.Position = emp(k).Col(indexwCol).Position;
                        GlobalBest.Velocity = emp(k).Col(indexwCol).Velocity;
                        GlobalBest.Cost = emp(k).Col(indexwCol).Cost;
                        GlobalBest.Power = emp(k).Col(indexwCol).Power;
                        emp(k).Col(indexwCol).Position = temp.Position;
                        emp(k).Col(indexwCol).Velocity = temp.Velocity;
                        emp(k).Col(indexwCol).Cost = temp.Cost;
                        emp(k).Col(indexwCol).Power = temp.Power;
                    end

                    % Exchange a colony with its local imperialist if the colony has more power than its local imperialist
                    if emp(k).Col(indexwCol).Power > emp(k).Imp.Power
                        [emp(k).Imp, emp(k).Col(indexwCol)] = deal(emp(k).Col(indexwCol), emp(k).Imp);                
                    end
            end
        end
    end
        BestPower(it) = GlobalBest.Power; 
     
    if FEs >= MaxFEs
        break;
    end
    end
    FEsRuns(nRun) = FEs;
    % Save Global Bests
    CostRuns(nRun) = GlobalBest.Cost;
end

BestResults = min(CostRuns);
MeanResults = mean(CostRuns);
MedianResults = median(CostRuns);
SDResults = std(CostRuns);

%% Show results
disp(['To reduce compute time on f1(Sphere), the MaxFEs was set to ' num2str(MaxFEs)]);
disp('With this adjustment, the results are :');
disp([ ' Best = '  num2str(BestResults)]);
disp([ ' Mean = '  num2str(MeanResults)]);
disp([ ' Median = '  num2str(MedianResults)]);
disp([ ' SD = '  num2str(SDResults)]);

%% Save results
save('BestResultsf1.mat','BestResults');
save('MeanResultsf1.mat','MeanResults');
save('MedianResultsf1.mat','MedianResults');
save('SDResultsf1.mat','SDResults');

%% For convergence speed analysis, we calculated the total FEs considering the stopping criterion of "global best's power tolerance<=10^-2"
%BestFEsRunha = min(FEsRuns);
%MeanFEsRunha = mean(FEsRuns);
toc