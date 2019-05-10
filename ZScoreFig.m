%% Load Data
Directory = '/Users/barabasi/Dropbox/Shared_Folders/Brain-Determistic-Dani/Code/RandAnalyzed/';
ChemAnalyzedName = 'CeChemRandBCs';
GapAnalyzedName = 'NewCGapBCs';
FullAnalyzedName = 'CeFullRandBCs';
EndCap = '.txt';

howManyReps = 1000;

RandChem_SourceSizes = zeros(91,howManyReps);
RandChem_DestSizes = zeros(63,howManyReps);
RandChem_SourcePlusDest = zeros(91+63,howManyReps);
RandChem_SourceTimesDest = zeros(200,howManyReps);
RandChem_SizesMat = zeros(91,63,howManyReps);
RandChemNumBCsNeur = zeros(2,279,howManyReps);
RandChemPercBCsNeur = zeros(2,279,howManyReps);

RandGap_SourceSizes = zeros(91,howManyReps);
RandGap_DestSizes = zeros(63,howManyReps);
RandGap_SourcePlusDest = zeros(91+63,howManyReps);
RandGap_SourceTimesDest = zeros(200,howManyReps);
RandGap_SizesMat = zeros(91,63,howManyReps);
RandGapNumBCsNeur = zeros(2,279,howManyReps);
RandGapPercBCsNeur = zeros(2,279,howManyReps);

% RandFull_SourceSizes = zeros(91,howManyReps);
% RandFull_DestSizes = zeros(63,howManyReps);
% RandFull_SourcePlusDest = zeros(91+63,howManyReps);
% RandFull_SourceTimesDest = zeros(200,howManyReps);
% RandFull_SizesMat = zeros(91,63,howManyReps);

BioGap_SourceSizes = zeros(91,1);
BioGap_DestSizes = zeros(63,1);
BioGap_SourcePlusDest = zeros(91+63,1);
BioGap_SourceTimesDest = zeros(200,1);
BioGap_SizesMat = zeros(91,63);
BioGapNumBCsNeur = zeros(2,279,howManyReps);
BioGapPercBCsNeur = zeros(2,279,howManyReps);

BioChem_SourceSizes = zeros(91,1);
BioChem_DestSizes = zeros(63,1);
BioChem_SourcePlusDest = zeros(91+63,1);
BioChem_SourceTimesDest = zeros(200,1);
BioChem_SizesMat = zeros(91,63);
BioChemNumBCsNeur = zeros(2,279,howManyReps);
BioChemPercBCsNeur = zeros(2,279,howManyReps);

% BioFull_SourceSizes = zeros(91,1);
% BioFull_DestSizes = zeros(63,1);
% BioFull_SourcePlusDest = zeros(91+63,1);
% BioFull_SourceTimesDest = zeros(200,1);
% BioFull_SizesMat = zeros(91,63);

% Number of BCs per Bio Chem Node
% A = fopen('/Users/barabasi/Dropbox/Shared_Folders/Brain-Determistic-Dani/Code/ChemBCImproved.txt');
% 
% out = '';
% while ~strcmp(out,'List of maximal bicliques:')
%     out = fgetl(A);
% end
% out = fgetl(A);
% out = fgetl(A);
% countNum =1;
% 
% while ~strcmp(out,'')
%     C = strsplit(out,'|');
%     uSize = length(strfind(C{1},' '));
%     vSize = length(strfind(C{2},' ')) -1;
%     BioChem_SizesMat(uSize,vSize) = BioChem_SizesMat(uSize,vSize)+1;
%     BioChem_SourceSizes(uSize) = BioChem_SourceSizes(uSize)+1;
%     BioChem_DestSizes(vSize) = BioChem_DestSizes(vSize) +1;
%     BioChem_SourcePlusDest(uSize+vSize) = BioChem_SourcePlusDest(uSize+vSize)+1;
%     BioChem_SourceTimesDest(uSize*vSize) =BioChem_SourceTimesDest(uSize*vSize) +1;
%         SourceNodesStore = regexp(C{1},'\d*','Match');
%     SourceNodes = zeros(1,length(SourceNodesStore));
%     for i = 1:length(SourceNodesStore)
%         SourceNodes(i) = str2double(SourceNodesStore{i});
%     end
%     BioChemNumBCsNeur(1,SourceNodes) = BioChemNumBCsNeur(1,SourceNodes)+1;
%     
%     vSize = length(strfind(C{2},' ')) -1;
%     DestNodesStore = regexp(C{2},'\d*','Match');
%     DestNodes = zeros(1,length(DestNodesStore));
%     for i = 1:length(DestNodesStore)
%         DestNodes(i) = str2double(DestNodesStore{i})-279;
%     end
%     BioChemNumBCsNeur(2,DestNodes) = BioChemNumBCsNeur(2,DestNodes)+1;
%     countNum = countNum+1;
%     out = fgetl(A);
% end
% BioChemPercBCsNeur(:,:) = BioChemNumBCsNeur(:,:)/countNum;
% 
% 
% 
% % Number of BCs per Bio Gap Node
% A = fopen('/Users/barabasi/Dropbox/Shared_Folders/Brain-Determistic-Dani/Code/CElegansOld/CEO_Gap_BCs.txt');
% out = '';
% while ~strcmp(out,'List of maximal bicliques:')
%     out = fgetl(A);
% end
% out = fgetl(A);
% out = fgetl(A);
% countNum =1;
% 
% while ~strcmp(out,'')
%     C = strsplit(out,'|');
%     uSize = length(strfind(C{1},' '));
%     vSize = length(strfind(C{2},' ')) -1;
%     BioGap_SizesMat(uSize,vSize) = BioGap_SizesMat(uSize,vSize)+1;
%     BioGap_SourceSizes(uSize) = BioGap_SourceSizes(uSize)+1;
%     BioGap_DestSizes(vSize) = BioGap_DestSizes(vSize) +1;
%     BioGap_SourcePlusDest(uSize+vSize) = BioGap_SourcePlusDest(uSize+vSize)+1;
%     BioGap_SourceTimesDest(uSize*vSize) =BioGap_SourceTimesDest(uSize*vSize) +1;
%     SourceNodesStore = regexp(C{1},'\d*','Match');
%     SourceNodes = zeros(1,length(SourceNodesStore));
%     for i = 1:length(SourceNodesStore)
%         SourceNodes(i) = str2double(SourceNodesStore{i});
%     end
%     BioGapNumBCsNeur(1,SourceNodes) = BioGapNumBCsNeur(1,SourceNodes)+1;
%     
%     DestNodesStore = regexp(C{2},'\d*','Match');
%     DestNodes = zeros(1,length(DestNodesStore));
%     for i = 1:length(DestNodesStore)
%         DestNodes(i) = str2double(DestNodesStore{i})-279;
%     end
%     BioGapNumBCsNeur(2,DestNodes) = BioGapNumBCsNeur(2,DestNodes)+1;
%     countNum = countNum+1;
%         out = fgetl(A);
% end
% BioGapPercBCsNeur(:,:) = BioGapNumBCsNeur(:,:)/countNum;

% % Number of BCs per Bio Full Node
% A = fopen('/Users/barabasi/Dropbox/Shared_Folders/Brain-Determistic-Dani/Code/LineageCElegans.txt');
% out = '';
% while ~strcmp(out,'List of maximal bicliques:')
%     out = fgetl(A);
% end
% out = fgetl(A);
% countNum =1;
% 
% while ~strcmp(out,'')
%     C = strsplit(out,'|');
%     uSize = length(strfind(C{1},' '));
%     vSize = length(strfind(C{2},' ')) -1;
%     BioFull_SizesMat(uSize,vSize) = BioFull_SizesMat(uSize,vSize)+1;
%     BioFull_SourceSizes(uSize) = BioFull_SourceSizes(uSize)+1;
%     BioFull_DestSizes(vSize) = BioFull_DestSizes(vSize) +1;
%     BioFull_SourcePlusDest(uSize+vSize) = BioFull_SourcePlusDest(uSize+vSize)+1;
%     BioFull_SourceTimesDest(uSize*vSize) =BioFull_SourceTimesDest(uSize*vSize) +1;
%     countNum = countNum+1;
%     out = fgetl(A);
% end



for repNum = 1:howManyReps
    % Number of BCs per Rand Chem Node
%     A = fopen(strcat(Directory,ChemAnalyzedName,num2str(repNum),EndCap));
%     out = '';
%     while ~strcmp(out,'List of maximal bicliques:')
%         out = fgetl(A);
%     end
%     out = fgetl(A);
%     out = fgetl(A);
%     
%     countNum =1;
%     
%     while ~strcmp(out,'')
%         C = strsplit(out,'|');
%         uSize = length(strfind(C{1},' '));
%         vSize = length(strfind(C{2},' ')) -1;
%         RandChem_SizesMat(uSize,vSize,repNum) = RandChem_SizesMat(uSize,vSize,repNum)+1;
%         RandChem_SourceSizes(uSize,repNum) = RandChem_SourceSizes(uSize,repNum)+1;
%         RandChem_DestSizes(vSize,repNum) = RandChem_DestSizes(vSize,repNum) +1;
%         RandChem_SourcePlusDest(uSize+vSize,repNum) = RandChem_SourcePlusDest(uSize+vSize,repNum)+1;
%         RandChem_SourceTimesDest(uSize*vSize,repNum) = RandChem_SourceTimesDest(uSize*vSize,repNum) +1;
%         SourceNodesStore = regexp(C{1},'\d*','Match');
%         SourceNodes = zeros(1,length(SourceNodesStore));
%         for i = 1:length(SourceNodesStore)
%             SourceNodes(i) = str2double(SourceNodesStore{i});
%         end
%         RandChemNumBCsNeur(1,SourceNodes,repNum) = RandChemNumBCsNeur(1,SourceNodes,repNum)+1;
%         
%         DestNodesStore = regexp(C{2},'\d*','Match');
%         DestNodes = zeros(1,length(DestNodesStore));
%         for i = 1:length(DestNodesStore)
%             DestNodes(i) = str2double(DestNodesStore{i})-279;
%         end
%         RandChemNumBCsNeur(2,DestNodes,repNum) = RandChemNumBCsNeur(2,DestNodes,repNum)+1;
%         countNum = countNum+1;
%         out = fgetl(A);
%     end
%     RandChemPercBCsNeur(:,:,repNum) = RandChemNumBCsNeur(:,:,repNum)/countNum;

    %Number of BCs per Rand Gap Node
    A = fopen(strcat(Directory,GapAnalyzedName,num2str(repNum),EndCap));
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
        vSize = length(strfind(C{2},' ')) -1;
        RandGap_SizesMat(uSize,vSize,repNum) = RandGap_SizesMat(uSize,vSize,repNum)+1;
        RandGap_SourceSizes(uSize,repNum) = RandGap_SourceSizes(uSize,repNum)+1;
        RandGap_DestSizes(vSize,repNum) = RandGap_DestSizes(vSize,repNum) +1;
        RandGap_SourcePlusDest(uSize+vSize,repNum) = RandGap_SourcePlusDest(uSize+vSize,repNum)+1;
        RandGap_SourceTimesDest(uSize*vSize,repNum) =RandGap_SourceTimesDest(uSize*vSize,repNum) +1;
        SourceNodesStore = regexp(C{1},'\d*','Match');
        SourceNodes = zeros(1,length(SourceNodesStore));
        for i = 1:length(SourceNodesStore)
            SourceNodes(i) = str2double(SourceNodesStore{i});
        end
        RandGapNumBCsNeur(1,SourceNodes,repNum) = RandGapNumBCsNeur(1,SourceNodes,repNum)+1;
        
        vSize = length(strfind(C{2},' ')) -1;
        DestNodesStore = regexp(C{2},'\d*','Match');
        DestNodes = zeros(1,length(DestNodesStore));
        for i = 1:length(DestNodesStore)
            DestNodes(i) = str2double(DestNodesStore{i})-279;
        end
        RandGapNumBCsNeur(2,DestNodes,repNum) = RandGapNumBCsNeur(2,DestNodes,repNum)+1;
        countNum = countNum+1;
        out = fgetl(A);
    end
        RandGapPercBCsNeur(:,:,repNum) = RandGapNumBCsNeur(:,:,repNum)/countNum;

    
%     %Number of BCs per Rand Full Node
%     A = fopen(strcat(Directory,FullAnalyzedName,num2str(repNum),EndCap));
%     out = '';
%     while ~strcmp(out,'List of maximal bicliques:')
%         out = fgetl(A);
%     end
%     out = fgetl(A);
%     countNum =1;
%     
%     while ~strcmp(out,'')
%         C = strsplit(out,'|');
%         uSize = length(strfind(C{1},' '));
%         vSize = length(strfind(C{2},' ')) -1;
%         RandFull_SizesMat(uSize,vSize,repNum) = RandFull_SizesMat(uSize,vSize,repNum)+1;
%         RandFull_SourceSizes(uSize,repNum) = RandFull_SourceSizes(uSize,repNum)+1;
%         RandFull_DestSizes(vSize,repNum) = RandFull_DestSizes(vSize,repNum) +1;
%         RandFull_SourcePlusDest(uSize+vSize,repNum) = RandFull_SourcePlusDest(uSize+vSize,repNum)+1;
%         RandFull_SourceTimesDest(uSize*vSize,repNum) =RandFull_SourceTimesDest(uSize*vSize,repNum) +1;
%         countNum = countNum+1;
%         out = fgetl(A);
%     end
    
end

%%
ZCols = Zs;
ZCols(Zs>2) = 2;
ZCols(Zs<-2) = 0;
ZCols(and(Zs>-2,Zs<2)) = -1;
ZCols(isnan(Zs)) = 1;

%%
imagesc(1.5:5.5,1.5:5.5,ZCols(2:6,2:6));colormap(flipud(flag))
xticks(1.5:5.5)
xticklabels(2:6)
yticks(1.5:5.5)
yticklabels(2:6)
for i = 2:5
    line([0,6],[i,i])
    line([i,i],[0,6])
end

%%
minVal = 2;
maxVal = 50;
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