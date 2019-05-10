%%
T1 = readtable('/Users/barabasi/Dropbox/Shared_Folders/Brain-Determistic-Dani/Code/NeuronLineage_Part1.xls');
T2 = readtable('/Users/barabasi/Dropbox/Shared_Folders/Brain-Determistic-Dani/Code/NeuronLineage_Part2.xls');
A = table2cell(T1);
B = table2cell(T2);
LineageIdent = unique([A(:,1);A(:,2);B(:,1);B(:,2)]);
clear T1 T2 A B;

T = readtable('/Users/barabasi/Dropbox/Shared_Folders/Brain-Determistic-Dani/Code/CElegansOld/NeuronConnect.xls');
T.Nbr = [];
T = T(cellfun(@isempty, strfind(T.Type, 'NMJ')), :);

TChem = T(cellfun(@isempty, strfind(T.Type, 'EJ')), :);
TChem = TChem(cellfun(@isempty, strfind(TChem.Type, 'R')), :);
TChem.Type = [];

cEdgeList = strrep(table2cell(TChem),' ',''); %Convert from table and remove white space
clear TChem;
CellInd = unique(cEdgeList);
LineageChemMat = zeros(length(LineageIdent),length(LineageIdent));
subSize = size(LineageChemMat);
for i = 1:length(cEdgeList)
    source = find(strcmp(LineageIdent, cEdgeList{i,1}));
    dest = find(strcmp(LineageIdent, cEdgeList{i,2}));
    if and(~isempty(source),~isempty(dest))
        LineageChemMat(sub2ind(subSize,source,dest)) = 1;
    end
end
clear source dest i subSize cEdgeList A B

TGap = T(cellfun(@isempty, strfind(T.Type, 'R')), :);
TGap = TGap(cellfun(@isempty, strfind(TGap.Type, 'S')), :);
TGap.Type = [];

cEdgeList = strrep(table2cell(TGap),' ',''); %Convert from table and remove white space
clear TGap;
LineageGapMat = zeros(length(LineageIdent),length(LineageIdent));
subSize = size(LineageGapMat);
for i = 1:length(cEdgeList)
    source = find(strcmp(LineageIdent, cEdgeList{i,1}));
    dest = find(strcmp(LineageIdent, cEdgeList{i,2}));
    if and(~isempty(source),~isempty(dest))
        LineageGapMat(sub2ind(subSize,source,dest)) = 1;
    end
end
clear source dest i subSize cEdgeList A B
LineageFullMat = double((LineageGapMat + LineageChemMat)>0);

save('CEO_Connectomes','LineageFullMat','LineageGapMat','LineageChemMat');
%% Save Connectomes and Rands
% Graph = double(LineageChemMat);
% Name = '/Users/barabasi/Dropbox/Shared_Folders/Brain-Determistic-Dani/Code/CElegansOld/CEO_Chem.txt';
% dlmwrite(Name,full(Graph),' ')
% 
% Graph = double(LineageGapMat);
% Name = '/Users/barabasi/Dropbox/Shared_Folders/Brain-Determistic-Dani/Code/CElegansOld/CEO_Gap.txt';
% dlmwrite(Name,full(Graph),' ')
% 
% Graph = double(LineageFullMat);
% Name = '/Users/barabasi/Dropbox/Shared_Folders/Brain-Determistic-Dani/Code/CElegansOld/CEO_Full.txt';
% dlmwrite(Name,full(Graph),' ')

for i = 1:1000
    rng(i)
%     Graph1 = dir_generate_srand(LineageChemMat);
%     Name = strcat(strcat('/Users/barabasi/Dropbox/Shared_Folders/Brain-Determistic-Dani/Code/CElegansOld/CElegansOldRands/CEO_Chem_Rand',int2str(i)),'.txt');
%     dlmwrite(Name,full(Graph1),' ')
%     
    Graph2 = sym_generate_srand(LineageGapMat);
    Name = strcat(strcat('/Users/barabasi/Dropbox/Shared_Folders/Brain-Determistic-Dani/Code/CEN_Gap/CEN_Gap_Rand',int2str(i)),'.txt');
    dlmwrite(Name,full(Graph2),' ')
    
%     Graph = double((Graph1 + Graph2)>0);
%     Name = strcat(strcat('/Users/barabasi/Dropbox/Shared_Folders/Brain-Determistic-Dani/Code/CElegansOld/CElegansOldRands/CEO_Full_Rand',int2str(i)),'.txt');
%     dlmwrite(Name,full(Graph),' ')
end