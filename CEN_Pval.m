%% Indexing:
% 1 and 2: Gene Source and Dest
% 3 and 4: Lineage Source and Dest
% 5 and 6: Class Source and Dest
% 7 and 8: Neurotransmitter
% 9 and 10: Type

%% Node IDs
T1 = readtable('/Users/barabasi/Dropbox/Shared_Folders/Brain-Determistic-Dani/Code/NeuronLineage_Part1.xls');
T2 = readtable('/Users/barabasi/Dropbox/Shared_Folders/Brain-Determistic-Dani/Code/NeuronLineage_Part2.xls');
A = table2cell(T1);
B = table2cell(T2);
LineageDists = [A;B];
LineageIdent = unique([A(:,1);A(:,2);B(:,1);B(:,2)]);
clear T1 T2 A B;

%% Type
T = table2cell(readtable('/Users/barabasi/Dropbox/Shared_Folders/Brain-Determistic-Dani/Code/node_type.csv'));
NameTypeCorrespondance = T(ismember(T(:,1),LineageIdent),:);
[~, ix] = sort(NameTypeCorrespondance(:,1));
NameTypeCorrespondance = NameTypeCorrespondance(ix,:);

%% Class
T = readtable('CE_Class_NTM.xlsx');
T.Notes = [];
T.Neurotransmitter = [];
NameClassCorrespondance = table2cell(T);
NameClassCorrespondance = NameClassCorrespondance(ismember(NameClassCorrespondance(:,2),LineageIdent),:);
[~,~,ix] = unique(NameClassCorrespondance(:,1));
[numOfClass,b] = hist(ix,unique(ix));
[HowManyOfClassSize,ClassSize] = hist(numOfClass,unique(numOfClass));

%% Neurotransmitter
T = readtable('/Users/barabasi/Dropbox/Shared_Folders/Brain-Determistic-Dani/Code/CE_Class_NTM.xlsx');
T.Notes = [];
T.NeuronClass = [];
NameTransCorrespondance = table2cell(T);
NameTransCorrespondance = NameTransCorrespondance(ismember(NameTransCorrespondance(:,1),LineageIdent),:);

%% Gene Expression
load('/Users/barabasi/Dropbox/Shared_Folders/Brain-Determistic-Dani/Code/GeneExpressionNew.mat');
allPossiblePerms = nchoosek(1:279,2);
thetas = zeros(1,length(allPossiblePerms));

for i = 1:length(allPossiblePerms)
    neur1 = allPossiblePerms(i,1);
    neur2 = allPossiblePerms(i,2);
    neur1Exp = expression(neur1,:);
    neur2Exp = expression(neur2,:);
    n_11 = length(intersect(find(neur1Exp),find(neur2Exp)));
    n_00 = length(intersect(find(neur1Exp==0),find(neur2Exp==0)));
    n_01 = length(intersect(find(neur1Exp==0),find(neur2Exp==1)));
    n_10 = length(intersect(find(neur1Exp==1),find(neur2Exp==0)));
    
    sums = [n_10+n_11,n_00+n_01,n_10+n_00,n_11+n_01];
    sums(find(sums == 0)) = 1;
    
    norm = sqrt(prod(sums));
    
    r_theta = (n_11*n_00-n_01*n_10)/norm;
    thetas(i) = r_theta;
end

%% Lineage Dists
% Create Distance Matrix for Lineages
LineageDistMat = zeros(279,279);
for i = 1:length(LineageDists)
    LineageDistMat(find(strcmp(LineageIdent, LineageDists{i,1})),find(strcmp(LineageIdent, LineageDists{i,2}))) = LineageDists{i,3};
end

% Get Lineage Distance Distribution
dists = cell2mat(LineageDists(:,3));

meanDists = mean(dists);
stdDists = std(dists);

%%
A = fopen('/Users/barabasi/Dropbox/Shared_Folders/Brain-Determistic-Dani/Code/ChemBCImproved.txt');

out = '';
while ~strcmp(out,'Number of maximal bicliques:')
    out = fgetl(A);
end
out = fgetl(A);
pStorer = nan(10,str2double(out));

while ~strcmp(out,'List of maximal bicliques:')
    out = fgetl(A);
end
out = fgetl(A);
out = fgetl(A);

countNum =1;
warning('off')
while ~strcmp(out,'')
    C = strsplit(out,'|');
    
    %% Source
    uSize = length(strfind(C{1},' '));
    SourceNodes = str2double(regexp(C{1},'\d*','Match'));
    
    %Type
    SourceClass = NameTypeCorrespondance(SourceNodes,2);
    [test1,~,ub] = unique(SourceClass);
    SourceCounts = histc(ub, 1:length(test1));
    [maxCount,maxInd] = max(SourceCounts);
    numMaxClass = sum(strcmp(NameTypeCorrespondance(:,2), test1{maxInd}));
    pStorer(9,countNum) = exp(sum(log(1:(279-numMaxClass)))+sum(log(1:(279-uSize)))+sum(log(1:uSize))+sum(log(1:numMaxClass))-sum(log(1:(279-numMaxClass-uSize+maxCount)))-sum(log(1:279))-sum(log(1:maxCount))-sum(log(1:(numMaxClass-maxCount)))-sum(log(1:(uSize-maxCount))));
    
    % Class
    SourceClass = NameClassCorrespondance(SourceNodes,1);
    [test1,~,ub] = unique(SourceClass);
    SourceCounts = histc(ub, 1:length(test1));
    [maxCount,maxInd] = max(SourceCounts);
    if maxCount > 1
        numMaxClass = sum(strcmp(NameClassCorrespondance(:,1), test1{maxInd}));
        probNumerator = exp(sum(log(1:(279-numMaxClass)))+sum(log(1:(279-uSize)))+sum(log(1:uSize))+sum(log(1:numMaxClass))-sum(log(1:(279-numMaxClass-uSize+maxCount)))-sum(log(1:279))-sum(log(1:maxCount))-sum(log(1:(numMaxClass-maxCount)))-sum(log(1:(uSize-maxCount))));
        pStorer(5,countNum) = HowManyOfClassSize(find(ClassSize == numMaxClass))*probNumerator;
    end
    
    % Neurotransmitter  
    SourceClass = NameTransCorrespondance(SourceNodes,2);
    [test1,~,ub] = unique(SourceClass);
    SourceCounts = histc(ub, 1:length(test1));
    [maxCount,maxInd] = max(SourceCounts);
    numMaxClass = sum(strcmp(NameTransCorrespondance(:,2), test1{maxInd}));
    pStorer(7,countNum) = exp(sum(log(1:(279-numMaxClass)))+sum(log(1:(279-uSize)))+sum(log(1:uSize))+sum(log(1:numMaxClass))-sum(log(1:(279-numMaxClass-uSize+maxCount)))-sum(log(1:279))-sum(log(1:maxCount))-sum(log(1:(numMaxClass-maxCount)))-sum(log(1:(uSize-maxCount))));
    
    %% Destination
    vSize = length(strfind(C{2},' ')) -1;
    DestNodes = str2double(regexp(C{2},'\d*','Match'));
    
    
    % Type
    DestClass = NameTypeCorrespondance(DestNodes-279,2);
    [test2,~,ub] = unique(DestClass);
    DestCounts = histc(ub, 1:length(test2));
    [maxCount,maxInd] = max(DestCounts);
    numMaxClass = sum(strcmp(NameTypeCorrespondance(:,2), test2{maxInd}));
    pStorer(10,countNum) = exp(sum(log(1:(279-numMaxClass)))+sum(log(1:(279-vSize)))+sum(log(1:vSize))+sum(log(1:numMaxClass))-sum(log(1:(279-numMaxClass-vSize+maxCount)))-sum(log(1:279))-sum(log(1:maxCount))-sum(log(1:(numMaxClass-maxCount)))-sum(log(1:(vSize-maxCount))));
    
    % Class
    DestClass = NameClassCorrespondance(DestNodes-279,1);
    [test2,~,ub] = unique(DestClass);
    DestCounts = histc(ub, 1:length(test2));
    [maxCount,maxInd] = max(DestCounts);
    if maxCount > 1
        numMaxClass = sum(strcmp(NameClassCorrespondance(:,1), test2{maxInd}));
        probNumerator = exp(sum(log(1:(279-numMaxClass)))+sum(log(1:(279-vSize)))+sum(log(1:vSize))+sum(log(1:numMaxClass))-sum(log(1:(279-numMaxClass-vSize+maxCount)))-sum(log(1:279))-sum(log(1:maxCount))-sum(log(1:(numMaxClass-maxCount)))-sum(log(1:(vSize-maxCount))));
        pStorer(6,countNum) = HowManyOfClassSize(find(ClassSize == numMaxClass))*probNumerator;
    end
    
    % Neurotransmitter
    DestClass = NameTransCorrespondance(DestNodes-279,2);
    [test2,~,ub] = unique(DestClass);
    DestCounts = histc(ub, 1:length(test2));
    [maxCount,maxInd] = max(DestCounts);
    numMaxClass = sum(strcmp(NameTransCorrespondance(:,2), test2{maxInd}));
    pStorer(8,countNum) = exp(sum(log(1:(279-numMaxClass)))+sum(log(1:(279-vSize)))+sum(log(1:vSize))+sum(log(1:numMaxClass))-sum(log(1:(279-numMaxClass-vSize+maxCount)))-sum(log(1:279))-sum(log(1:maxCount))-sum(log(1:(numMaxClass-maxCount)))-sum(log(1:(vSize-maxCount))));
    
    % Genes
    if and(uSize>1,vSize>1)
        NodeCombos = nchoosek(SourceNodes,2);
        rThetaBtwnSourceNodes = zeros(1,size(NodeCombos,1));
        for comboLoop = 1:size(NodeCombos,1)
            
            neur1 = NodeCombos(comboLoop,1);
            neur2 = NodeCombos(comboLoop,2);
            neur1Exp = expression(neur1,:);
            neur2Exp = expression(neur2,:);
            n_11 = length(intersect(find(neur1Exp),find(neur2Exp)));
            n_00 = length(intersect(find(neur1Exp==0),find(neur2Exp==0)));
            n_01 = length(intersect(find(neur1Exp==0),find(neur2Exp==1)));
            n_10 = length(intersect(find(neur1Exp==1),find(neur2Exp==0)));
            
            sums = [n_10+n_11,n_00+n_01,n_10+n_00,n_11+n_01];
            sums(find(sums == 0)) = 1;
            
            norm = sqrt(prod(sums));
            
            r_theta = (n_11*n_00-n_01*n_10)/norm;
            
            rThetaBtwnSourceNodes(1,comboLoop) = r_theta;
        end
        [~,pStorer(1,countNum),~] = kstest2(rThetaBtwnSourceNodes,thetas,'tail','smaller');
        
        NodeCombos = nchoosek(DestNodes - 279,2);
        rthetaBtwnDestNodes = zeros(1,size(NodeCombos,1));
        for comboLoop = 1:size(NodeCombos,1)
            neur1 = NodeCombos(comboLoop,1);
            neur2 = NodeCombos(comboLoop,2);
            neur1Exp = expression(neur1,:);
            neur2Exp = expression(neur2,:);
            n_11 = length(intersect(find(neur1Exp),find(neur2Exp)));
            n_00 = length(intersect(find(neur1Exp==0),find(neur2Exp==0)));
            n_01 = length(intersect(find(neur1Exp==0),find(neur2Exp==1)));
            n_10 = length(intersect(find(neur1Exp==1),find(neur2Exp==0)));
            
            sums = [n_10+n_11,n_00+n_01,n_10+n_00,n_11+n_01];
            sums(find(sums == 0)) = 1;
            
            norm = sqrt(prod(sums));
            
            r_theta = (n_11*n_00-n_01*n_10)/norm;
            
            rthetaBtwnDestNodes(1,comboLoop) = r_theta;
        end
        [~,pStorer(2,countNum),~] = kstest2(rthetaBtwnDestNodes,thetas,'tail','smaller');
    end
    
    % Lineage
    if and(uSize>1,vSize>1)
        NodeCombos = nchoosek(SourceNodes,2);
        distBtwnSourceNodes = zeros(1,size(NodeCombos,1));
        for comboLoop = 1:size(NodeCombos,1)
            distBtwnSourceNodes(1,comboLoop) = LineageDistMat(NodeCombos(comboLoop,1),NodeCombos(comboLoop,2));
        end
        [~,pStorer(3,countNum),~] = kstest2(distBtwnSourceNodes,dists,'tail','larger');
        
        NodeCombos = nchoosek(DestNodes - 279,2);
        distBtwnDestNodes = zeros(1,size(NodeCombos,1));
        for comboLoop = 1:size(NodeCombos,1)
            distBtwnDestNodes(1,comboLoop) = LineageDistMat(NodeCombos(comboLoop,1),NodeCombos(comboLoop,2));
        end
        [~,pStorer(4,countNum),~] = kstest2(distBtwnDestNodes,dists,'tail','larger');
        
    end
    
    countNum = countNum+1;
    out = fgetl(A);
end
warning('on')