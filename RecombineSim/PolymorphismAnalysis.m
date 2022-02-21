clear
vars = readtable('SK1_SNP_WeightMap','Delimiter','\t');
vars = table2array(vars);
vars(:,3) = 1;
fid1 = fopen('ChrSizesS288cH4L2_L2HG.txt','r');
chrsizes = round(cell2mat(textscan(fid1, '%d','HeaderLines',1))/100);
data = readtable('SK1_VarAnalysisInput','Delimiter','\t');
uID = table2array(unique(data(:,1)));
for s=1:length(uID)
    indices = find(strcmp(data{:,{'Genotype'}},uID(s)));
    chr_data = data(indices,:);
    chr_data = table2array(chr_data(:,2:3));
    for j=1:16
        varmap = zeros(1,chrsizes(j));
        polym = find(vars(:,1)==j);
        chr_var = vars(polym,:);
        chr_var(:,2) = round(chr_var(:,2)/100,0);
        idx = find(chr_data(:,1)==j);
        pos = chr_data(idx,:);
        pos(:,2) = round(pos(:,2)/100,0);
        bin=20;
        pos(pos(:,2)+bin>length(varmap),:) = [];
        for wg=1:length(chr_var)
            varmap(chr_var(wg,2)) = varmap(chr_var(wg,2))+chr_var(wg,3);
        end
        dens=0;
        for co=1:length(pos)
            dens(co) = sum(varmap((pos(co,2)-bin):(pos(co,2)+bin)));
        end
        res{s,j}=dens;
    end
end
hold on
for k=1:length(uID)
    ksdensity(cell2mat(res(k,:)),[0:100],'Function','cdf','Bandwidth',1); %Plot
    [x,y] = ksdensity(cell2mat(res(k,:)),[0:100],'Function','cdf','Bandwidth',1); %Data
    x=x';
    x(:,2)=y;
    cdf{k}=x;
end
xlim([0,20])
res=cell2mat(cdf);
res=array2table(res);
s=0;
for b=1:2:size(res,2)
    s=s+1;
    res.Properties.VariableNames{b} = uID{s};
end
writetable(res,"Results.txt",'Delimiter','\t')