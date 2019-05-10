pCutoff = 0.05;
sourceSize = 8;
destSize = 4;

%% Import AdjMat
T1 = readtable('/Users/barabasi/Dropbox/Shared_Folders/Brain-Determistic-Dani/Code/NeuronLineage_Part1.xls');
T2 = readtable('/Users/barabasi/Dropbox/Shared_Folders/Brain-Determistic-Dani/Code/NeuronLineage_Part2.xls');
A = table2cell(T1);
B = table2cell(T2);
LineageIdent = unique([A(:,1);A(:,2);B(:,1);B(:,2)]);
TNull = readtable('/Users/barabasi/Dropbox/Shared_Folders/Brain-Determistic-Dani/Code/herm_full_edgelist.csv');
T = TNull(cellfun(@isempty, strfind(TNull.Type, 'electrical')), :);
T.Weight = [];
T.Type = [];
cEdgeList = strrep(table2cell(T),' ',''); %Convert from table and remove white space
CellInd = unique(cEdgeList);

LineageFullMat = sparse(length(LineageIdent),length(LineageIdent));
LineageChemMat = sparse(length(LineageIdent),length(LineageIdent));
LineageGapMat = sparse(length(LineageIdent),length(LineageIdent));

subSize = size(LineageFullMat);
for i = 1:length(cEdgeList)
    source = find(strcmp(LineageIdent, cEdgeList{i,1}));
    dest = find(strcmp(LineageIdent, cEdgeList{i,2}));
    if and(~isempty(source),~isempty(dest))
        LineageFullMat(sub2ind(subSize,source,dest)) = 1;
        LineageChemMat(sub2ind(subSize,source,dest)) = 1;
        
    end
end


T = TNull(cellfun(@isempty, strfind(TNull.Type, 'chemical')), :);
T.Weight = [];
T.Type = [];
cEdgeList = strrep(table2cell(T),' ',''); %Convert from table and remove white space
for i = 1:length(cEdgeList)
    source = find(strcmp(LineageIdent, cEdgeList{i,1}));
    dest = find(strcmp(LineageIdent, cEdgeList{i,2}));
    if and(~isempty(source),~isempty(dest))
        LineageFullMat(sub2ind(subSize,source,dest)) = 1;
        LineageFullMat(sub2ind(subSize,dest,source)) = 1;
        
        LineageGapMat(sub2ind(subSize,source,dest)) = 1;
        LineageGapMat(sub2ind(subSize,dest,source)) = 1;
    end
end

LineageMat = LineageChemMat;
clear T T1 T2 TNull cEdgeList A B

T = readtable('CE_Class_NTM.xlsx');
T.Notes = [];
T.Neurotransmitter = [];
NameClassCorrespondance = table2cell(T);
NameClassCorrespondance = NameClassCorrespondance(ismember(NameClassCorrespondance(:,2),LineageIdent),:);

T = table2cell(readtable('node_type.csv'));
NameTypeCorrespondance = T(ismember(T(:,1),LineageIdent),:);
[~, ix] = sort(NameTypeCorrespondance(:,1));
NameTypeCorrespondance = NameTypeCorrespondance(ix,:);
[~,~,Type_ix] = unique(NameTypeCorrespondance(:,2));
Type_Cols = [1 0 0 ; 0 1 0; 0 0 1];

T = readtable('CE_Class_NTM.xlsx');
T.Notes = [];
T.NeuronClass = [];
NameTransCorrespondance = table2cell(T);
NameTransCorrespondance = NameTransCorrespondance(ismember(NameTransCorrespondance(:,1),LineageIdent),:);
[~,~,Trans_ix] = unique(NameTransCorrespondance(:,2));
Trans_Cols = linspecer(max(Trans_ix));

%% Get BiClique Sizes and Distribution


A = fopen('ChemBCImproved.txt');

out = '';
while ~strcmp(out,'Number of maximal bicliques:')
    out = fgetl(A);
end
out = fgetl(A);
SizeStorer = zeros(2,str2double(out));

while ~strcmp(out,'List of maximal bicliques:')
    out = fgetl(A);
end
out = fgetl(A);
out = fgetl(A);
countNum =1;

BCmat2 = zeros(91,63);

while ~strcmp(out,'')
    C = strsplit(out,'|');
    uSize = length(strfind(C{1},' '));
    vSize = length(strfind(C{2},' ')) -1;
    SizeStorer(:,countNum) = [uSize vSize];
    BCmat2(uSize,vSize) = BCmat2(uSize,vSize)+1;
    countNum = countNum+1;
    out = fgetl(A);
end
clear out uSize vSize A C
%%
load('ChemPValues.mat')
iDs = find(and(SizeStorer(1,:) == sourceSize,SizeStorer(2,:) == destSize));
pS = Chem_Pmeta(:,iDs);
pSRelevant = find(and(pS(3,:)<pCutoff,pS(4,:)<pCutoff));
whichBCs = iDs(pSRelevant);

%%
A = fopen('ChemBCImproved.txt');

numCheck =1;
BCPlotNum = whichBCs(numCheck);
out = '';
while ~strcmp(out,'List of maximal bicliques:')
    out = fgetl(A);
end
out = fgetl(A);
out = fgetl(A);
BCCount = 1;

while ~strcmp(out,'')
    
    if BCCount == BCPlotNum
        C = strsplit(out,'|');
        
        SourceNodesStore = regexp(C{1},'\d*','Match');
        SourceNodes = zeros(1,length(SourceNodesStore));
        for i = 1:length(SourceNodesStore)
            SourceNodes(i) = str2double(SourceNodesStore{i});
        end
        
        DestNodesStore = regexp(C{2},'\d*','Match');
        DestNodes = zeros(1,length(DestNodesStore));
        for i = 1:length(DestNodesStore)
            DestNodes(i) = str2double(DestNodesStore{i});
        end
        hF = subplot(1,2,1);
        hold off
        H = plot(digraph(BipartizeSubGraph(LineageMat(SourceNodes,DestNodes-length(LineageMat)))),'Layout','layered','NodeLabel',[],'ShowArrows','off');
        hold on
        title(strcat(num2str(sourceSize),'x',num2str(destSize),' Biclique'));
        classMems = NameClassCorrespondance([SourceNodes,DestNodes-279],1);
        TypeMems = Type_ix([SourceNodes,DestNodes-279],1);
        TransMems = Trans_ix([SourceNodes,DestNodes-279],1);
        [~,~,icOls] = unique(NameClassCorrespondance([SourceNodes,DestNodes-279],1));
        cols = linspecer(max(icOls));
        for i=1:(length(SourceNodes)+length(DestNodes))
            if i >length(SourceNodes)
                text(H.XData(i),H.YData(i)-0.03,LineageIdent{DestNodes(i-length(SourceNodes))-279},'fontsize',13,'Rotation',-90);
                
            else
                text(H.XData(i),H.YData(i)+0.03,LineageIdent{SourceNodes(i)},'fontsize',13,'Rotation',90);
            end
            highlight(H,i,'NodeColor',cols(icOls(i),:),'MarkerSize',12)
        end
        
        sourcePs={'Gene Expression: ',num2str(pS(3,pSRelevant(numCheck)));'Lineage Distance: ',num2str(pS(1,pSRelevant(numCheck)));'Class Enrichment: ',num2str(pS(5,pSRelevant(numCheck)));'Location Enrichment: ',num2str(pS(9,pSRelevant(numCheck)));'Neurotransmitter Enrichment: ',num2str(pS(11,pSRelevant(numCheck)));'Type Enrichment: ',num2str(pS(7,pSRelevant(numCheck)))};
        destPs={'Gene Expression: ',num2str(pS(4,pSRelevant(numCheck)));'Lineage Distance: ',num2str(pS(2,pSRelevant(numCheck)));'Class Enrichment: ',num2str(pS(6,pSRelevant(numCheck)));'Location Enrichment: ',num2str(pS(10,pSRelevant(numCheck)));'Neurotransmitter Enrichment: ',num2str(pS(12,pSRelevant(numCheck)));'Type Enrichment: ',num2str(pS(8,pSRelevant(numCheck)))};
        set(hF,'color',[1 1 1])
        hA=axes; set(hA,'color',[1 1 1],'visible','off')
        %text(0.5,0.9,sourcePs(:,1),'FontSize',15,'FontWeight','bold')
        text(0.5,1+0.05,'Source','FontSize',16,'FontWeight','bold','Color','b')
        text(0.5,0.23+0.05,'Destination','FontSize',16,'FontWeight','bold','Color','b')
        
        for i = 1:6
            if str2double(sourcePs(i,2))<0.05
                text(0.9,1-(i-1)*0.054,sourcePs(i,2),'FontSize',15,'FontWeight','bold','Color','r')
                text(0.5,1-(i-1)*0.054,sourcePs(i,1),'FontSize',15,'FontWeight','bold','Color','r')
                
            else
                text(0.9,1-(i-1)*0.054,sourcePs(i,2),'FontSize',15,'FontWeight','bold')
                text(0.5,1-(i-1)*0.054,sourcePs(i,1),'FontSize',15,'FontWeight','bold')
                
            end
            if str2double(destPs(i,2))<0.05
                text(0.9,0.23-(i-1)*0.054,destPs(i,2),'FontSize',15,'FontWeight','bold','Color','r')
                text(0.5,0.23-(i-1)*0.054,destPs(i,1),'FontSize',15,'FontWeight','bold','Color','r')
            else
                text(0.9,0.23-(i-1)*0.054,destPs(i,2),'FontSize',15,'FontWeight','bold')
                text(0.5,0.23-(i-1)*0.054,destPs(i,1),'FontSize',15,'FontWeight','bold')
            end
        end
        %text(0.5,0.1,destPs(:,1),'FontSize',15,'FontWeight','bold')
        k = waitforbuttonpress;
        numCheck = numCheck +1;
        BCPlotNum = whichBCs(numCheck);
    end
    BCCount = BCCount + 1;
    out = fgetl(A);
end

%% Gap
T1 = readtable('/Users/barabasi/Dropbox/Shared_Folders/Brain-Determistic-Dani/Code/NeuronLineage_Part1.xls');
T2 = readtable('/Users/barabasi/Dropbox/Shared_Folders/Brain-Determistic-Dani/Code/NeuronLineage_Part2.xls');
A = table2cell(T1);
B = table2cell(T2);
LineageIdent = unique([A(:,1);A(:,2);B(:,1);B(:,2)]);
TNull = readtable('/Users/barabasi/Dropbox/Shared_Folders/Brain-Determistic-Dani/Code/herm_full_edgelist.csv');
T = TNull(cellfun(@isempty, strfind(TNull.Type, 'electrical')), :);
T.Weight = [];
T.Type = [];
cEdgeList = strrep(table2cell(T),' ',''); %Convert from table and remove white space
CellInd = unique(cEdgeList);

LineageFullMat = sparse(length(LineageIdent),length(LineageIdent));
LineageChemMat = sparse(length(LineageIdent),length(LineageIdent));
LineageGapMat = sparse(length(LineageIdent),length(LineageIdent));

subSize = size(LineageFullMat);
for i = 1:length(cEdgeList)
    source = find(strcmp(LineageIdent, cEdgeList{i,1}));
    dest = find(strcmp(LineageIdent, cEdgeList{i,2}));
    if and(~isempty(source),~isempty(dest))
        LineageFullMat(sub2ind(subSize,source,dest)) = 1;
        LineageChemMat(sub2ind(subSize,source,dest)) = 1;
        
    end
end


T = TNull(cellfun(@isempty, strfind(TNull.Type, 'chemical')), :);
T.Weight = [];
T.Type = [];
cEdgeList = strrep(table2cell(T),' ',''); %Convert from table and remove white space
for i = 1:length(cEdgeList)
    source = find(strcmp(LineageIdent, cEdgeList{i,1}));
    dest = find(strcmp(LineageIdent, cEdgeList{i,2}));
    if and(~isempty(source),~isempty(dest))
        LineageFullMat(sub2ind(subSize,source,dest)) = 1;
        LineageFullMat(sub2ind(subSize,dest,source)) = 1;
        
        LineageGapMat(sub2ind(subSize,source,dest)) = 1;
        LineageGapMat(sub2ind(subSize,dest,source)) = 1;
    end
end

LineageMat = LineageGapMat;
clear T T1 T2 TNull cEdgeList A B

T = readtable('CE_Class_NTM.xlsx');
T.Notes = [];
T.Neurotransmitter = [];
NameClassCorrespondance = table2cell(T);
NameClassCorrespondance = NameClassCorrespondance(ismember(NameClassCorrespondance(:,2),LineageIdent),:);

T = table2cell(readtable('node_type.csv'));
NameTypeCorrespondance = T(ismember(T(:,1),LineageIdent),:);
[~, ix] = sort(NameTypeCorrespondance(:,1));
NameTypeCorrespondance = NameTypeCorrespondance(ix,:);
[~,~,Type_ix] = unique(NameTypeCorrespondance(:,2));
Type_Cols = [1 0 0 ; 0 1 0; 0 0 1];

T = readtable('CE_Class_NTM.xlsx');
T.Notes = [];
T.NeuronClass = [];
NameTransCorrespondance = table2cell(T);
NameTransCorrespondance = NameTransCorrespondance(ismember(NameTransCorrespondance(:,1),LineageIdent),:);
[~,~,Trans_ix] = unique(NameTransCorrespondance(:,2));
Trans_Cols = linspecer(max(Trans_ix));



A = fopen('GapBCImproved.txt');

out = '';
while ~strcmp(out,'Number of maximal bicliques:')
    out = fgetl(A);
end
out = fgetl(A);
SizeStorer = zeros(2,str2double(out));

while ~strcmp(out,'List of maximal bicliques:')
    out = fgetl(A);
end
out = fgetl(A);
out = fgetl(A);
countNum =1;

BCmat2 = zeros(91,63);

while ~strcmp(out,'')
    C = strsplit(out,'|');
    uSize = length(strfind(C{1},' '));
    vSize = length(strfind(C{2},' ')) -1;
    SizeStorer(:,countNum) = [uSize vSize];
    BCmat2(uSize,vSize) = BCmat2(uSize,vSize)+1;
    countNum = countNum+1;
    out = fgetl(A);
end
clear out uSize vSize A C

load('GapPValues.mat')
iDs = find(and(SizeStorer(1,:) == sourceSize,SizeStorer(2,:) == destSize));
pS = Gap_Pmeta(:,iDs);
pSRelevant = find(and(pS(3,:)<pCutoff,pS(4,:)<pCutoff));
whichBCs = iDs(pSRelevant);

A = fopen('GapBCImproved.txt');

numCheck =1;
BCPlotNum = whichBCs(numCheck);
out = '';
while ~strcmp(out,'List of maximal bicliques:')
    out = fgetl(A);
end
out = fgetl(A);
out = fgetl(A);
BCCount = 1;

while ~strcmp(out,'')
    
    if BCCount == BCPlotNum
        C = strsplit(out,'|');
        
        SourceNodesStore = regexp(C{1},'\d*','Match');
        SourceNodes = zeros(1,length(SourceNodesStore));
        for i = 1:length(SourceNodesStore)
            SourceNodes(i) = str2double(SourceNodesStore{i});
        end
        
        DestNodesStore = regexp(C{2},'\d*','Match');
        DestNodes = zeros(1,length(DestNodesStore));
        for i = 1:length(DestNodesStore)
            DestNodes(i) = str2double(DestNodesStore{i});
        end
        hF = subplot(1,2,1);
        hold off
        H = plot(digraph(BipartizeSubGraph(LineageMat(SourceNodes,DestNodes-length(LineageMat)))),'Layout','layered','NodeLabel',[],'ShowArrows','off');
        hold on
        title(strcat(num2str(sourceSize),'x',num2str(destSize),' Biclique'));
        classMems = NameClassCorrespondance([SourceNodes,DestNodes-279],1);
        TypeMems = Type_ix([SourceNodes,DestNodes-279],1);
        TransMems = Trans_ix([SourceNodes,DestNodes-279],1);
        [~,~,icOls] = unique(NameClassCorrespondance([SourceNodes,DestNodes-279],1));
        cols = linspecer(max(icOls));
        for i=1:(length(SourceNodes)+length(DestNodes))
            if i >length(SourceNodes)
                text(H.XData(i),H.YData(i)-0.03,LineageIdent{DestNodes(i-length(SourceNodes))-279},'fontsize',13,'Rotation',-90);
                
            else
                text(H.XData(i),H.YData(i)+0.03,LineageIdent{SourceNodes(i)},'fontsize',13,'Rotation',90);
            end
            highlight(H,i,'NodeColor',cols(icOls(i),:),'MarkerSize',12)
        end
        
        sourcePs={'Gene Expression: ',num2str(pS(3,pSRelevant(numCheck)));'Lineage Distance: ',num2str(pS(1,pSRelevant(numCheck)));'Class Enrichment: ',num2str(pS(5,pSRelevant(numCheck)));'Location Enrichment: ',num2str(pS(9,pSRelevant(numCheck)));'Neurotransmitter Enrichment: ',num2str(pS(11,pSRelevant(numCheck)));'Type Enrichment: ',num2str(pS(7,pSRelevant(numCheck)))};
        destPs={'Gene Expression: ',num2str(pS(4,pSRelevant(numCheck)));'Lineage Distance: ',num2str(pS(2,pSRelevant(numCheck)));'Class Enrichment: ',num2str(pS(6,pSRelevant(numCheck)));'Location Enrichment: ',num2str(pS(10,pSRelevant(numCheck)));'Neurotransmitter Enrichment: ',num2str(pS(12,pSRelevant(numCheck)));'Type Enrichment: ',num2str(pS(8,pSRelevant(numCheck)))};
        set(hF,'color',[1 1 1])
        hA=axes; set(hA,'color',[1 1 1],'visible','off')
        %text(0.5,0.9,sourcePs(:,1),'FontSize',15,'FontWeight','bold')
        text(0.5,1+0.05,'Source','FontSize',16,'FontWeight','bold','Color','b')
        text(0.5,0.23+0.05,'Destination','FontSize',16,'FontWeight','bold','Color','b')
        
        for i = 1:6
            if str2double(sourcePs(i,2))<0.05
                text(0.9,1-(i-1)*0.054,sourcePs(i,2),'FontSize',15,'FontWeight','bold','Color','r')
                text(0.5,1-(i-1)*0.054,sourcePs(i,1),'FontSize',15,'FontWeight','bold','Color','r')
                
            else
                text(0.9,1-(i-1)*0.054,sourcePs(i,2),'FontSize',15,'FontWeight','bold')
                text(0.5,1-(i-1)*0.054,sourcePs(i,1),'FontSize',15,'FontWeight','bold')
                
            end
            if str2double(destPs(i,2))<0.05
                text(0.9,0.23-(i-1)*0.054,destPs(i,2),'FontSize',15,'FontWeight','bold','Color','r')
                text(0.5,0.23-(i-1)*0.054,destPs(i,1),'FontSize',15,'FontWeight','bold','Color','r')
            else
                text(0.9,0.23-(i-1)*0.054,destPs(i,2),'FontSize',15,'FontWeight','bold')
                text(0.5,0.23-(i-1)*0.054,destPs(i,1),'FontSize',15,'FontWeight','bold')
            end
        end
        %text(0.5,0.1,destPs(:,1),'FontSize',15,'FontWeight','bold')
        k = waitforbuttonpress;
        numCheck = numCheck +1;
        BCPlotNum = whichBCs(numCheck);
    end
    BCCount = BCCount + 1;
    out = fgetl(A);
end