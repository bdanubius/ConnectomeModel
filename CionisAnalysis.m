%% Get Mats
CAdjRaw = double(xlsread('AdjMatRaw.xlsx')>0);
[~,~,SourceNames] = xlsread('RawNames.xlsx','A1:A205');
SourceNames=cellfun(@num2str,SourceNames,'un',0);
[~,~,DestNames] = xlsread('RawNames.xlsx','B1:B215');
DestNames=cellfun(@num2str,DestNames,'un',0);
[RemainingNames,SourceI,DestI] = intersect(SourceNames,DestNames);
CSquare = CAdjRaw(SourceI,DestI);
%[~,j] = find(CSquare(find(sum(CSquare,2)==1),:));
%eliminate = [find(CSquare(:,intersect(j,find(sum(CSquare,2)==0)))),find(and(sum(CSquare,1).'==0,sum(CSquare,2)==0)),intersect(j,find(sum(CSquare,2)==0))];
eliminate = find(and(sum(CSquare,1).'==0,sum(CSquare,2)==0));
CSquare(eliminate,:) = [];
CSquare(:,eliminate) = [];
RemainingNames(eliminate) = [];
%% Save New
Graph = double(CSquare);
Name = '/Users/barabasi/Dropbox/Shared_Folders/Brain-Determistic-Dani/Code/CionisIntestinalis/CionIntestMat.txt';
dlmwrite(Name,full(Graph),' ')

for i = 1:1000
    rng(i)
    Graph = dir_generate_srand(CSquare);
    Name = strcat(strcat('/Users/barabasi/Dropbox/Shared_Folders/Brain-Determistic-Dani/Code/CionisIntestinalis/CionisRands/CionisRand',int2str(i)),'.txt');
    dlmwrite(Name,full(Graph),' ')
end

%% Size of BCs
howManyReps = 1000;

RandChem_SizesMat = zeros(50,50,howManyReps);
RandNumBCsNeur = zeros(2,197,howManyReps);
RandPercBCsNeur = zeros(2,197,howManyReps);


BioChem_SizesMat = zeros(50,50);
BioNumBCsNeur = zeros(2,197);
BioPercBCsNeur = zeros(2,197);


% Number of BCs per Bio Chem Node
A = fopen('/Users/barabasi/Dropbox/Shared_Folders/Brain-Determistic-Dani/Code/CionisIntestinalis/CionBCs.txt');

out = '';
while ~strcmp(out,'List of maximal bicliques:')
    out = fgetl(A);
end
out = fgetl(A);
out = fgetl(A);
countNum =1;

while ~strcmp(out,'')
    C = strsplit(out,'|');
    uSize = length(strfind(C{1},' '));
    SourceNodesStore = regexp(C{1},'\d*','Match');
    SourceNodes = zeros(1,length(SourceNodesStore));
    for i = 1:length(SourceNodesStore)
        SourceNodes(i) = str2double(SourceNodesStore{i});
    end
    BioNumBCsNeur(1,SourceNodes) = BioNumBCsNeur(1,SourceNodes)+1;
    
    vSize = length(strfind(C{2},' ')) -1;
    DestNodesStore = regexp(C{2},'\d*','Match');
    DestNodes = zeros(1,length(DestNodesStore));
    for i = 1:length(DestNodesStore)
        DestNodes(i) = str2double(DestNodesStore{i})-197;
    end
    BioNumBCsNeur(2,DestNodes) = BioNumBCsNeur(2,DestNodes)+1;
    countNum = countNum+1;
    BioChem_SizesMat(uSize,vSize) = BioChem_SizesMat(uSize,vSize)+1;
    
    
    out = fgetl(A);
end
BioPercBCsNeur(:,:) = BioNumBCsNeur(:,:)/countNum;


for repNum = 1:howManyReps
    % Number of BCs per Rand Chem Node
    A = fopen(strcat('/Users/barabasi/Dropbox/Shared_Folders/Brain-Determistic-Dani/Code/CionisIntestinalis/CionisRands/CionRandBC',num2str(repNum),'.txt'));
    out = '';
    while ~strcmp(out,'List of maximal bicliques:')
        out = fgetl(A);
    end
    out = fgetl(A);
    out = fgetl(A);
    
    countNum =1;
    
    while ~strcmp(out,'')
        C = strsplit(out,'|');
        uSize = length(strfind(C{1},' '));
        SourceNodesStore = regexp(C{1},'\d*','Match');
        SourceNodes = zeros(1,length(SourceNodesStore));
        for i = 1:length(SourceNodesStore)
            SourceNodes(i) = str2double(SourceNodesStore{i});
        end
        RandNumBCsNeur(1,SourceNodes,repNum) = RandNumBCsNeur(1,SourceNodes,repNum)+1;
        
        vSize = length(strfind(C{2},' ')) -1;
        DestNodesStore = regexp(C{2},'\d*','Match');
        DestNodes = zeros(1,length(DestNodesStore));
        for i = 1:length(DestNodesStore)
            DestNodes(i) = str2double(DestNodesStore{i})-197;
        end
        RandNumBCsNeur(2,DestNodes,repNum) = RandNumBCsNeur(2,DestNodes,repNum)+1;
        RandChem_SizesMat(uSize,vSize,repNum) = RandChem_SizesMat(uSize,vSize,repNum)+1;
        countNum = countNum+1;
        out = fgetl(A);
    end
    RandPercBCsNeur(:,:,repNum) = RandNumBCsNeur(:,:,repNum)/countNum;
end

clearvars -except RandChem_SizesMat BioChem_SizesMat RandNumBCsNeur BioNumBCsNeur BioPercBCsNeur RandPercBCsNeur

%% Z-Score Fig
%save('Ciona1000RepsData')
Zs = (BioChem_SizesMat-mean(RandChem_SizesMat,3))./std(RandChem_SizesMat,0,3);
ZCols = Zs;
ZCols(Zs>2) = 2;
ZCols(Zs<-2) = 0;
ZCols(and(Zs>-2,Zs<2)) = -1;
ZCols(isnan(Zs)) = 1;

minVal = 2;
maxVal = 25;
imagesc((minVal-0.5):(maxVal-0.5),(minVal-0.5):(maxVal-0.5),ZCols(minVal:maxVal,minVal:maxVal)); colormap(flipud(flag))
xticks((minVal-0.5):(maxVal-0.5))
xticklabels(minVal:maxVal)
yticks((minVal-0.5):(maxVal-0.5))
yticklabels(minVal:maxVal)
set(gca,'xaxisLocation','top')

for i = (minVal-1):maxVal
    line([0,maxVal],[i,i],'Color','w')
    line([i,i],[0,maxVal],'Color','w')
end
pbaspect([1 1 1])

%% Degree Vs BCs
figure
subplot(1,2,1)
scatter(sum(CSquare,2),BioPercBCsNeur(1,:))
xlim([0,50])
xlabel('k')
ylabel('Percent of BCs')
hold on
scatter(sum(CSquare,2),mean(RandPercBCsNeur(1,:,:),3))
subplot(1,2,2)
scatter(sum(CSquare,1),BioPercBCsNeur(2,:))
xlabel('k')
ylabel('Percent of BCs')
hold on
scatter(sum(CSquare,1),mean(RandPercBCsNeur(2,:,:),3))
xlim([0,50])