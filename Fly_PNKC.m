load('Fly_PN_KC.mat','unnamed');
Fly_PNKCAll = unnamed;

Graph = double(Fly_PNKCAll>4);
Name = '/Users/barabasi/Dropbox/Shared_Folders/Brain-Determistic-Dani/Code/Flies/PNKC.txt';
dlmwrite(Name,full(Graph),' ')
ROP = Fly_PNKCAll>4;
%% Save Files
for i = 1:1000
    rng(i)
    
    Graph = dir_generate_srand(ROP);
    Name = strcat(strcat('/Users/barabasi/Dropbox/Shared_Folders/Brain-Determistic-Dani/Code/Flies/PNKC/PNKCRand',int2str(i)),'.txt');
    dlmwrite(Name,full(Graph),' ')
end

%% Size of BCs
%howManyReps = 1000;

% RandChem_SizesMat = zeros(96,96,howManyReps);
% 
% 
% BioChem_SizesMat = zeros(96,96);
% 
% 
% % Number of BCs per Bio Chem Node
% A = fopen('/Users/barabasi/Dropbox/Shared_Folders/Brain-Determistic-Dani/Code/Flies/PNKCBCs.txt');
% 
% out = '';
% while ~strcmp(out,'List of maximal bicliques:')
%     out = fgetl(A);
% end
% out = fgetl(A);
% out = fgetl(A);
% 
% countNum =1;
% 
% while ~strcmp(out,'')
%     C = strsplit(out,'|');
%     uSize = length(strfind(C{1},' '));
%     SourceNodesStore = regexp(C{1},'\d*','Match');
%     SourceNodes = zeros(1,length(SourceNodesStore));
%     for i = 1:length(SourceNodesStore)
%         SourceNodes(i) = str2double(SourceNodesStore{i});
%     end
%     
%     vSize = length(strfind(C{2},' ')) -1;
%     DestNodesStore = regexp(C{2},'\d*','Match');
%     DestNodes = zeros(1,length(DestNodesStore));
%     for i = 1:length(DestNodesStore)
%         DestNodes(i) = str2double(DestNodesStore{i})-96;
%     end
%     countNum = countNum+1;
%     BioChem_SizesMat(uSize,vSize) = BioChem_SizesMat(uSize,vSize)+1;
%     
%     
%     out = fgetl(A);
% end


for repNum = 944:howManyReps
    % Number of BCs per Rand Chem Node
    A = fopen(strcat('/Users/barabasi/Dropbox/Shared_Folders/Brain-Determistic-Dani/Code/Flies/PNKC/PNKCRandBC',num2str(repNum),'.txt'));
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
        
        vSize = length(strfind(C{2},' ')) -1;
        DestNodesStore = regexp(C{2},'\d*','Match');
        DestNodes = zeros(1,length(DestNodesStore));
        for i = 1:length(DestNodesStore)
            DestNodes(i) = str2double(DestNodesStore{i})-96;
        end
        RandChem_SizesMat(uSize,vSize,repNum) = RandChem_SizesMat(uSize,vSize,repNum)+1;
        countNum = countNum+1;
        out = fgetl(A);
        
    end
    repNum
end

clearvars -except RandChem_SizesMat BioChem_SizesMat RandNumBCsNeur BioNumBCsNeur BioPercBCsNeur RandPercBCsNeur

%%

Zs = (BioChem_SizesMat-mean(RandChem_SizesMat,3))./std(RandChem_SizesMat,0,3);
ZCols = Zs;
ZCols(Zs>2) = 2;
ZCols(Zs<-2) = 0;
ZCols(and(Zs>-2,Zs<2)) = -1;
ZCols(isnan(Zs)) = 1;

minVal = 2;
maxVal = 47;
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