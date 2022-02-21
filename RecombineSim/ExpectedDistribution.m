clear
vars = readtable('SK1_SNP_WeightMap','Delimiter','\t');
vars = table2array(vars);
vars(:,3) = 1;
fid1 = fopen('ChrSizesS288cH4L2_L2HG.txt','r');
chrsizes = round(cell2mat(textscan(fid1, '%d','HeaderLines',1))/100);
data = readtable('DSB_WeightMap','Delimiter','\t');
data = table2array(data);
for j=1:16
    varmap = zeros(1,chrsizes(j));
    polym = find(vars(:,1)==j);
    chr_var = vars(polym,:);
    chr_var(:,2) = round(chr_var(:,2)/100,0);
    idx = find(data(:,1)==j);
    weights = data(idx,:);
    weights(:,2) = round(weights(:,2)/100,0);
    bin=5;
    weights(weights(:,2)+bin>length(varmap),:) = [];
    weights(find(weights(:,2)<=bin+1),:) = [];
    dsbmap = zeros(1,chrsizes(j));
    for wg=1:length(chr_var)
        varmap(chr_var(wg,2)) = varmap(chr_var(wg,2))+chr_var(wg,3);
    end
    for wg=1:length(weights)
        dsbmap(weights(wg,2)) = dsbmap(weights(wg,2))+weights(wg,3);
    end
    pos = randsample([1:length(dsbmap)], 1000, true,dsbmap);
    dens=0;
    for co=1:length(pos)
        dens(co) = sum(varmap((pos(co)-bin):(pos(co)+bin)));
    end
    res{j}=dens;
end
res=cell2mat(res);
ksdensity(res,[0:100],'Function','cdf','Bandwidth',1); %Plot
[x,y] = ksdensity(res,[0:100],'Function','cdf','Bandwidth',1); %Data
x=x';
x(:,2)=y;
dlmwrite("SK1_Expected_500bp.txt",x,'Delimiter','\t')
