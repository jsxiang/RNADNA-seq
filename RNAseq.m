classdef RNAseq < handle
    % to do
    % collapse mutants - done, in python
    % handle replicates - done
    % graph structure
    % find lengths
    % plot cleavage & identify switches - done
    % volcano plot - done
    % display information for specified sequences - done
    % need a way to calculate decay constant instead
    
    
  properties
    counts;
    seqs;
    refs;
    ctrls;
    val;
    quantctls;
    mutantfamily;

 end

  methods
    function obj=RNAseq()
    end

    function load(data,file)
        % Backward compatibility with old data.mat files
        if nargin<2
            file='combinedcounts.csv';
        end
        fprintf('Loading data from %s...',file);
        d=csvread(file,0,1);
        fid=fopen(file);

        s=textscan(fid,'%s');
        fclose(fid);

        seqs={};
        for j=1:length(s{1})
            o=regexp(s{1}{j},'(\W*)([A-Z]+[-]?[A-Z]+)(\W*)','tokens');
            try
            seqs{end+1}=o{1}{2};
            catch
            seqs{end+1}='naw';

            end
        end
            data.counts=d;
            data.seqs=seqs;

        fprintf('done\n');
%         fn=setdiff(fieldnames(data),'pts');
%         for i=1:length(fn)
%             data.(fn{i})=data.(fn{i});
%         end
        %       obj.cleanup();
    end
    
    function loadfromprevious(data,file)
        fprintf('Loading data from %s....',file);
        d=load(file);
        fnload=fieldnames(d.data);
        for i=1:length(fnload)
           data.(fnload{i})=d.data.(fnload{i});
        end
    end
    
    function findctrls(data,ctrls)
        if nargin<2
            ctrls.seqs={'GCTGTCACCGGATGTGCTTTCCGGTCTGATGAGTCCGTGAGGACGAAACAGC',
                        'GCTGTCACCGGATGTGCTTTCCGGTACGTGAGGTCCGTGAGGACGAAACAGC',
                        'GCTGTCACCGGATGTTTCCGGTCTGATGAGTCCATAAGGACGAAACAGC',
                        'GCTGTCACCGGATACTTCCGGTCTGATGAGTCCCAGAGGACGAAACAGC',
                        'GCTGTCACCGGATGCATCCGGTCTGATGAGTCCCGCGGGACGAAACAGC'};
            ctrls.names={'sTRSV','sTRSVctl','g814','g833','g862'};
        end
        data.ctrls=ctrls;
        ctrlids=[];
        for i=1:length(data.ctrls.seqs)
            id=find(~cellfun('isempty',regexp(data.seqs,strcat('^',data.ctrls.seqs{i},'$'))));
            if ~isempty(id)
                ctrlids(i)=id;
            end
        end
        data.ctrls.ids=ctrlids;
    end
    
    function findvalids(data,val) % validated theo switches
        if nargin<2
            val.seqs={
                    'GCTGTCACCGGAATACCAGCATCGTCTTGATGCCCTTGGAAGTCCGGTCTGATGAGTCCCATAAGGACGAAACAGC',
                    'GCTGTCACCGGAATACCAGCATCGTCTTGATGCCCTTGGAAGTCCGGTCTGATGAGTCCAAAAAGGACGAAACAGC',
                    'GCTGTCACCGGAATACCAGCATCGTCTTGATGCCCTTGGAAGTCCGGTCTGATGAGTCCAGAAAGGACGAAACAGC',
                    'GCTGTCACCGGAATACCAGCATCGTCTTGATGCCCTTGGAAGTCCGGTCTGATGAGTCCCAGAAGGACGAAACAGC',
                    'GCTGTCACCGGAATACCAGCATCGTCTTGATGCCCTTGGAAGTCCGGTCTGATGAGTCCAGGAAGGACGAAACAGC',
                    'GCTGTCACCGGAATACCAGCATCGTCTTGATGCCCTTGGAAGTCCGGTCTGATGAGTCCAAAGAGGACGAAACAGC',
                    'GCTGTCACCGGAATACCAGCATCGTCTTGATGCCCTTGGAAGTCCGGTCTGATGAGTCCATTCAGGACGAAACAGC',
                    'GCTGTCACCGGAATACCAGCATCGTCTTGATGCCCTTGGAAGTCCGGTCTGATGAGTCCCTGAAGGACGAAACAGC',
                    'GCTGTCACCGGAATACCAGCATCGTCTTGATGCCCTTGGAAGTCCGGTCTGATGAGTCCCGGCAGGACGAAACAGC',
                    'GCTGTCACCGGAATACCAGCATCGTCTTGATGCCCTTGGAAGTCCGGTCTGATGAGTCCCTTGAGGACGAAACAGC',
                    'GCTGTCACCGGAATACCAGCATCGTCTTGATGCCCTTGGAAGTCCGGTCTGATGAGTCCTGGTAGGACGAAACAGC',
                    'GCTGTCACCGGAATACCAGCATCGTCTTGATGCCCTTGGAAGTCCGGTCTGATGAGTCCGAAAGGGACGAAACAGC',
                    };
             val.names={'CAUAA','AAAAA','AGAAA','CAGAA','AGGAA','AAAGA','AUUCA','CUGAA','CGGCA','CUUGA','UGGUA','GAAAG'};
%              val.names={'CAUAA','AAAAA','AGAAA','CAGAA','AGGAA','AAAGA','AUUCA','CGGCA','CUUGA','UGGUA','GAAAG'};
             
             
        end
        data.val=val;
        valids=[];
        for i=1:length(data.val.seqs)
            id=find(~cellfun('isempty',regexp(data.seqs,strcat('^',data.val.seqs{i},'$'))));
            if ~isempty(id)
                valids(i)=id;
            end
        end
        data.val.ids=valids;
    end
    
    function findquantctls(data,quantctls)
        if nargin<2
            quantctls.loop1={'GATGTGCTTTC',
                            'GAATACCAGCATCGTCTTGATGCCCTTGGAAGTC',
                            'GTGTATTACCCTAGTGGTCGAC',
                            'AGAAACTGTGAAGTATATCTTAAACCTGGGCACTTAAAAGATATATGGAGTTAGTAGTGCAACCTGCT',
                            'GTGCTTGGTACGTTATATTCAGC',
                            'GGTAGAACAGAATGCAAGTGCGTAGAAGCGCCTC'};
            quantctls.loop2={'ACAA',
                            'AGCT',
                            'GTGA',
                            'AGTT',
                            'CATA',
                            'ATTA',
                            'AGAA',
                            'TTAA',
                            'ACTA',
                            'TTGA'};

            quantctrlnames={'sTRSVctl','theoAAG','xan','cdGII','FA','zea'};
            rbz1='GCTGTCACCG';
            rbz2='CGGTCTGATGAGTCC';
            rbz3='GGACGAAACAGC';
            rbz3mut='GGACGAAACAGC';
            for i=1:length(quantctls.loop1)
                for j=1:length(quantctls.loop2)
                    quantctls.apt(i).seqs{j}=strcat(rbz1,quantctls.loop1{i},rbz2,quantctls.loop2{j},rbz3);
                    quantctls.apt(i).mutseqs{j}=strcat(rbz1,quantctls.loop1{i},rbz2,quantctls.loop2{j},rbz3mut);
                    quantctls.apt(i).name=quantctrlnames{i};
                end
            end
            quantctls.names=quantctrlnames;
        end
        data.quantctls=quantctls;
        for i=1:length(data.quantctls.apt)
            quantids=[];
            mutquantids=[];
            
            for j=1:length(data.quantctls.apt(i).seqs)
                id=find(~cellfun('isempty',regexp(data.seqs,strcat('^',data.quantctls.apt(i).seqs{j},'$'))));
                if ~isempty(id)
                    quantids(j)=id;
                end
                
                id=find(~cellfun('isempty',regexp(data.seqs,strcat('^',data.quantctls.apt(i).mutseqs{j},'$'))));
                if ~isempty(id)
                    mutquantids(j)=id;
                end 
            end
                data.quantctls.apt(i).ids=quantids;
                data.quantctls.apt(i).mutids=mutquantids;
        end
        
    end
    
    function findpairplots(data,ids,samplenames,figurename,aptaname)
    setfig(figurename);clf
    for i=1:length(ids)
        for j=1:length(ids)
        subplot(length(ids),length(ids),(i-1)*length(ids)+j)
        hold on
        x=data.counts(:,ids(i));
        y=data.counts(:,ids(j));
        xdata=x(x>2&y>2);
        ydata=y(x>2&y>2);
        plot(xdata,ydata,'.','MarkerSize',8)

        % plot controls

        for k=1:length(data.ctrls.ids)
            if data.ctrls.ids(k)>0
            plot(data.counts(data.ctrls.ids(k),ids(i)),data.counts(data.ctrls.ids(k),ids(j)),'+','MarkerSize',10)
            end
        end

        for k=1:length(data.val.ids)
            if data.val.ids(k)>0
            plot(data.counts(data.val.ids(k),ids(i)),data.counts(data.val.ids(k),ids(j)),'o','MarkerSize',10)
            end
        end
        
        if nargin>4
            cid=find(~cellfun('isempty',regexp(data.quantctls.names,aptaname)));
            plot(data.counts(data.quantctls.apt(cid).ids(data.quantctls.apt(cid).ids~=0),ids(i)),data.counts(data.quantctls.apt(cid).ids(data.quantctls.apt(cid).ids~=0),ids(j)),'d','MarkerSize',10)
            plot(data.counts(data.quantctls.apt(cid).mutids(data.quantctls.apt(cid).mutids~=0),ids(i)),data.counts(data.quantctls.apt(cid).mutids(data.quantctls.apt(cid).mutids~=0),ids(j)),'ks','MarkerSize',10)
        end

        x=linspace(min(data.counts(:,ids(i))),max(data.counts(:,ids(j))));
        plot(x,x,'k:','linewidth',1.5)
        xlabel(samplenames{ids(i)})
        ylabel(samplenames{ids(j)})
        set(gca,'fontsize',14)
        set(gca,'linewidth',1.5)
        set(gca,'YScale','log')
        set(gca,'XScale','log')

        axis([1 10^4.4 1 10^4.4])
        axis([2e1 1e4 2e1 1e4])
        grid off
        end

        if ((i-1)*length(ids)+j)==length(ids)*length(ids)
            
            legend('all seqs','sTRSV','sTRSVctl','g814','g833','g862',data.val.names{:},'1:1','location','best')
        end
    end
    end
    
    function plotreadcounts(data,i,j,figurename,outputfigname,fitrbz,filters)
        % i and j are scalars or vectors for the rounds/bc groups
        setfig(figurename);clf
        
        if nargin<7
            filters=1:length(data.counts(:,1));
        end
        scatter(sum(data.counts(filters,i),2),sum(data.counts(filters,j),2),'filled','MarkerFaceAlpha',0.3)
        [sum(data.counts(filters,i),2),sum(data.counts(filters,j),2)]
        % scatter(sum(RSII.counts(:,[10 11]),2),sum(RSII.counts(:,[1 2]),2),'filled','MarkerFaceAlpha',0.3)
        hold on

        for k=1:length(data.ctrls.ids)
            plot(sum(data.counts(data.ctrls.ids(k),i),2),sum(data.counts(data.ctrls.ids(k),j),2),'+','MarkerSize',14,'linewidth',2)
        %     plot(sum(RSII.counts(ctrlids(k),[10 11]),2),sum(RSII.counts(ctrlids(k),[1 2]),2),'+','MarkerSize',14,'linewidth',2)
        end

        if fitrbz
            mdl=fitlm(sum(data.counts(data.ctrls.ids,i),2),sum(data.counts(data.ctrls.ids,j),2));
            rsq=mdl.Rsquared.Adjusted;
        else
            mdl=fitlm(sum(data.counts(filters,i),2),sum(data.counts(filters,j),2));
            rsq=mdl.Rsquared.Adjusted;
        end

        fitcoeffs=[mdl.Coefficients.Estimate(1) mdl.Coefficients.Estimate(2)];
        x=linspace(0.2,max(sum(data.counts(:,j),2)));

        plot(x,fitcoeffs(1)+fitcoeffs(2)*x,'k:','linewidth',2)

        text(10,1000,sprintf('R^2=%0.3f',rsq),'fontsize',16)

    % plot(x,x,'k:','linewidth',2)

    %     title('sTRSVctl loopII:N4')
    %     xlabel('DNA read count rep1')
    %     ylabel('DNA read count rep2')

        set(gca,'fontsize',18)
        set(gca,'linewidth',2)
        set(gca,'YScale','log')
        set(gca,'XScale','log')
%         axis([1 10^3.4 1 10^3.4])
    %     legend('all seqs','sTRSV','sTRSVctl','g814','g833','g862','1:1','location','best')
        saveas(gca,outputfigname,'png')
    end
    
    function subrun=findRNADNAratio(data,RNAids1,RNAids2,DNAids1,DNAids2,figurename,aptamer,minthresh,plotcontrols,plotreplicates,plotgraph)
        normscalefactor=2;
        i=RNAids1;
        k=RNAids2;
        j=DNAids1;
        l=DNAids2;

        if nargin<4
            minthresh=100;
            plotcontrols='false';
            plotreplicates='true';
            plotgraph='true';
        end
        if plotgraph
            setfig(figurename)
        end
        % need to be clear that none of the counts are zero
        enoughreads=[data.counts(:,[RNAids1,RNAids2])>0 sum(data.counts(:,[DNAids1]),2)>minthresh sum(data.counts(:,[DNAids2]),2)>minthresh];
        
        hasaptamer=(~cellfun('isempty',regexp(data.seqs,aptamer)));
        goodstat=(sum(enoughreads,2)==4)&(hasaptamer');
        goodstat=(sum(enoughreads,2)==4);
        
        subrun.goodstat=goodstat;

        RNA_DNA1=log10(sum(data.counts(:,i),2)/sum(sum(data.counts(:,i),2))./sum(data.counts(:,j),2)*sum(sum(data.counts(:,j),2)));
        RNA_DNA2=log10(sum(data.counts(:,k),2)/sum(sum(data.counts(:,k),2))./sum(data.counts(:,l),2)*sum(sum(data.counts(:,l),2)));
        if ~plotcontrols & data.ctrls.ids(2)>0
            RNA_DNA1(data.ctrls.ids(2))=3;
            RNA_DNA2(data.ctrls.ids(2))=3;
        elseif ~plotcontrols
            data.ctrls.ids(2)=length(RNA_DNA1)+1;
            RNA_DNA1(data.ctrls.ids(2))=3;
            RNA_DNA2(data.ctrls.ids(2))=3;
            
        end
        RNA_DNA_good1=log10(sum(data.counts(goodstat,i),2)/sum(sum(data.counts(:,i),2))./sum(data.counts(goodstat,j),2)*sum(sum(data.counts(:,j),2)))-RNA_DNA1(data.ctrls.ids(2))+normscalefactor;
        RNA_DNA_good2=log10(sum(data.counts(goodstat,k),2)/sum(sum(data.counts(:,k),2))./sum(data.counts(goodstat,l),2)*sum(sum(data.counts(:,l),2)))-RNA_DNA2(data.ctrls.ids(2))+normscalefactor;

        
            
            
            
        subrun.seqs=data.seqs(goodstat);
        subrun.RDratio=[RNA_DNA_good1 RNA_DNA_good2];

        subrun.ctrls_RDratio=[RNA_DNA1(data.ctrls.ids)-RNA_DNA1(data.ctrls.ids(2))+normscalefactor RNA_DNA2(data.ctrls.ids)-RNA_DNA2(data.ctrls.ids(2))+normscalefactor];
        subrun.val_RDratio=[(RNA_DNA1(data.val.ids)-RNA_DNA1(data.ctrls.ids(2))+normscalefactor),(RNA_DNA2(data.val.ids)-RNA_DNA2(data.ctrls.ids(2))+normscalefactor)];
        ctrls_RDratio=[RNA_DNA1(data.ctrls.ids)-RNA_DNA1(data.ctrls.ids(2))+normscalefactor RNA_DNA2(data.ctrls.ids)-RNA_DNA2(data.ctrls.ids(2))+normscalefactor];
        if plotgraph

            scatter((RNA_DNA_good1),(RNA_DNA_good2),'filled','MarkerFaceAlpha',0.2,'MarkerFaceColor',[0.1 0.1 0.1]);
            scatter((RNA_DNA_good1),(RNA_DNA_good2),'filled','MarkerFaceAlpha',0.2);
            hold on

            for c=1:length(data.ctrls.ids)
                plot((RNA_DNA1(data.ctrls.ids(c))-RNA_DNA1(data.ctrls.ids(2))+normscalefactor),(RNA_DNA2(data.ctrls.ids(c))-RNA_DNA2(data.ctrls.ids(2))+normscalefactor),'+','MarkerSize',14,'linewidth',2)
                fprintf('%2.4f\t%2.4f\n',RNA_DNA1(data.ctrls.ids(c)),RNA_DNA2(data.ctrls.ids(c)))
            end

            for c=1:length(data.val.ids)
                plot((RNA_DNA1(data.val.ids(c))-RNA_DNA1(data.ctrls.ids(2))+normscalefactor),(RNA_DNA2(data.val.ids(c))-RNA_DNA2(data.ctrls.ids(2))+normscalefactor),'o','MarkerSize',14,'linewidth',2)
            end
        end

        x=linspace(min(RNA_DNA_good1),max(RNA_DNA_good1));
        % plot(x,x,':','linewidth',2)
        if plotcontrols
            % mdl=fitlm(sum(RS.motif(1).subrun(p).ctrls_RDratio,2),sum(RS.motif(1).subrun(p).ctrls_RDratio,2));
            mdl=fitlm((ctrls_RDratio(:,1)),ctrls_RDratio(:,2));
            % mdl=fitlm(RNA_DNA_good1,RNA_DNA_good2);
            rsq=mdl.Rsquared.Adjusted;

            fitcoeffs=[mdl.Coefficients.Estimate(1) mdl.Coefficients.Estimate(2)];
            if plotgraph
            plot(x,fitcoeffs(1)+fitcoeffs(2)*x,'k--','linewidth',2)
            end
        end

        if plotreplicates
            mdl=fitlm((RNA_DNA_good2),(RNA_DNA_good1));
            rsq=mdl.Rsquared.Adjusted;
% 
            fitcoeffs=[mdl.Coefficients.Estimate(1) mdl.Coefficients.Estimate(2)];
            if plotgraph
            plot(x,fitcoeffs(1)+fitcoeffs(2)*x,'--','linewidth',2)
            text(1.05,2.05,sprintf('R^2=%0.2f',rsq),'fontsize',16)
            end
            subrun.adjustedRDratio=[subrun.RDratio(:,1) fitcoeffs(1)+fitcoeffs(2)*subrun.RDratio(:,2)];
            subrun.adjustedRDratio_ctrls=[subrun.ctrls_RDratio(:,1) fitcoeffs(1)+fitcoeffs(2)*subrun.ctrls_RDratio(:,2)];
            if ~isempty(data.val.ids)
            subrun.adjustedRDratio_val=[subrun.val_RDratio(:,1) fitcoeffs(1)+fitcoeffs(2)*subrun.val_RDratio(:,2)];
            else
                subrun.adjustedRDratio_val=[];
            end
        end
        
      
        RNA_DNA_double=(sum(data.counts(:,[i k]),2)/sum(sum(data.counts(:,[i k]),2))./sum(data.counts(:,[j l]),2)*sum(sum(data.counts(:,[j l]),2)));
        RNA_DNA_double=(sum(data.counts(:,[i k]),2)./sum(data.counts(:,[j l]),2));
        
        subrun.goodstat=goodstat;
        subrun.totalcount=sum(data.counts(subrun.goodstat,[DNAids1 DNAids2]),2);
        subrun.DNAcount=[sum(data.counts(subrun.goodstat,[DNAids1 ]),2) sum(data.counts(subrun.goodstat,[ DNAids2]),2)];
        subrun.RNAcount=sum(data.counts(subrun.goodstat,[RNAids1 RNAids2]),2);
%         subrun.RDraw=sum(data.counts(subrun.goodstat,[RNAids1 RNAids2]),2)/sum(sum(data.counts(:,[RNAids1 RNAids2]),2))./sum(data.counts(subrun.goodstat,[DNAids1 DNAids2]),2)*sum(sum(data.counts(:,[DNAids1 DNAids2]),2));
        subrun.RDraw=RNA_DNA_double(subrun.goodstat);
        subrun.ctrls=data.ctrls;
        subrun.val=data.val;
        
        
        for l=1:length(data.ctrls.ids)
            try
                subrun.ctrls_RDraw(l)=RNA_DNA_double(data.ctrls.ids(l));

            catch
                subrun.ctrls_RDraw(l)=nan;
            end

            try
                subrun.ctrlcounts(l)=sum(data.counts(data.ctrls.ids(l),[DNAids1 DNAids2]),2);
            catch
                subrun.ctrlcounts(l)=nan;

            end
            
            
        end
        if ~isempty(data.val.ids)
        for l=1:length(data.val.ids)
            try
            subrun.valcounts(l)=sum(data.counts(data.val.ids(l),[DNAids1 DNAids2]),2);
            catch
            subrun.valcounts(l)=nan;
            end
        end
        end
%         xlabel(xnames{p})
%         ylabel(ynames{p})
        % xlabel(xnames)
        % ylabel(ynames)
        if plotgraph
        set(gca,'fontsize',20)
        set(gca,'linewidth',1.5)
        set(gca,'YScale','log')
        set(gca,'XScale','log')
        end

%         if p==2
%             legend('all seqs','sTRSV','sTRSVctl','g814','g833','g862','CATAA','AAAAA','AGAAA','CAGAA','AGGAA','linear fit','location','southeast')
%         end
% axis([0.8 3.5 0.8 3.5])

        
    end

    function plotsubrun(~,subrun,filters,adjusted)

        RNA_DNA_good1=subrun.RDratio(:,1);
        if nargin<4
            filters=1:length(RNA_DNA_good1);
            adjusted=true;
        end
        
        if adjusted
            RNA_DNA_good2=subrun.adjustedRDratio(:,2);
            ctrls_RDratio=subrun.adjustedRDratio_ctrls;
            val_RDratio=subrun.adjustedRDratio_val;
        else
            RNA_DNA_good2=subrun.RDratio(:,2);
            ctrls_RDratio=subrun.ctrls_RDratio;
            val_RDratio=subrun.val_RDratio;
        end

        scatter((RNA_DNA_good1(filters)),(RNA_DNA_good2(filters)),'filled','MarkerFaceAlpha',0.2);
        hold on

        for c=1:length(ctrls_RDratio(:,1))
            plot(ctrls_RDratio(c,1),ctrls_RDratio(c,2),'+','MarkerSize',14,'linewidth',2)
        end
        if ~isempty(subrun.val_RDratio)
        for c=1:length(subrun.val_RDratio(:,1))
            plot(val_RDratio(c,1),val_RDratio(c,2),'o','MarkerSize',14,'linewidth',2)
        end
        end

        x=linspace(min(RNA_DNA_good1),max(RNA_DNA_good1));
        % plot(x,x,':','linewidth',2)

        % mdl=fitlm(sum(RS.motif(1).subrun(p).ctrls_RDratio,2),sum(RS.motif(1).subrun(p).ctrls_RDratio,2));
%         mdl=fitlm(sum(subrun.ctrls_RDratio(:,1),2),sum(subrun.ctrls_RDratio(:,2),2));

        size(filters)
        size(RNA_DNA_good1)
        size(RNA_DNA_good2)
        mdl=fitlm(RNA_DNA_good1(filters),RNA_DNA_good2(filters));
        rsq=mdl.Rsquared.Adjusted;

        fitcoeffs=[mdl.Coefficients.Estimate(1) mdl.Coefficients.Estimate(2)];
        plot(x,fitcoeffs(1)+fitcoeffs(2)*x,'k--','linewidth',2)

        text(1,2.05,sprintf('R^2=%0.2f',rsq),'fontsize',16)
        text(1,1.85,sprintf('n=%d',length(filters)),'fontsize',16)

%         xlabel(xnames{p})
%         ylabel(ynames{p})
        % xlabel(xnames)
        % ylabel(ynames)
        set(gca,'fontsize',20)
        set(gca,'linewidth',1.5)
        set(gca,'YScale','log')
        set(gca,'XScale','log')

    end
    
    function combined=combinesubrun(~,subrun1,subrun2,filterbyaptamer,aptamer,alpha)
        
        [c,ia,ib]=intersect(subrun1.seqs,subrun2.seqs);
        
        RD1=subrun1.adjustedRDratio(ia,:);
        RD2=subrun2.adjustedRDratio(ib,:);
        
        ctrl1=subrun1.adjustedRDratio_ctrls;
        ctrl2=subrun2.adjustedRDratio_ctrls;
        
        val1=subrun1.adjustedRDratio_val;
        val2=subrun2.adjustedRDratio_val;

        RDcomb=[mean(RD1,2), mean(RD2,2)];
        ctrlcomb=[mean(ctrl1,2), mean(ctrl2,2)];
        valcomb=[mean(val1,2), mean(val2,2)];
        
        [Y1,I1]=sort(RDcomb(:,1));
        R1sorted=Y1(end:-5:(end-50));
        R2sorted=RDcomb(I1(end:-5:(end-50)),2);
        R1sorted=[];
        R2sorted=[];
        x1=[ctrl1(:,2);ctrl2(:,2)];
        x2=[ctrl1(:,1);ctrl2(:,1)];
        x1=[ctrlcomb(:,2); R2sorted];
        x2=[ctrlcomb(:,1); R1sorted];
        x1=[ctrlcomb([2 4 5],2); R2sorted]; % uncomment if cdG RNAseq
        x2=[ctrlcomb([2 4 5],1); R1sorted]; % uncomment if cdG RNAseq
        mdl=fitlm(x1,x2);
        fitcoeffs=[mdl.Coefficients.Estimate(1) mdl.Coefficients.Estimate(2)];

        adjustedRDcomb=[RDcomb(:,1) fitcoeffs(1)+fitcoeffs(2)*RDcomb(:,2)];
        adjustedctrlcomb=[ctrlcomb(:,1) fitcoeffs(1)+fitcoeffs(2)*ctrlcomb(:,2)];
        adjustedvalcomb=[valcomb(:,1) fitcoeffs(1)+fitcoeffs(2)*valcomb(:,2)];
                
%         adjustedRDcomb=RDcomb;
%         adjustedctrlcomb=ctrlcomb;
%         adjustedvalcomb=valcomb;
        
        
        if filterbyaptamer
        hasapt1=~cellfun('isempty',regexp(subrun1.seqs(ia),aptamer));
        else
            hasapt1=ones(length(adjustedRDcomb(:,1)),1);
        end
        
        scatter(adjustedRDcomb(hasapt1,1),adjustedRDcomb(hasapt1,2),70,'filled','markerfacealpha',0.2,'MarkerFaceColor',[0.1 0.1 0.1]);
        
        hold on
        
        combined.seqs=subrun1.seqs(ia(hasapt1));
        combined.RDratio=adjustedRDcomb(hasapt1,:);
        combined.ctrls_RDratio=adjustedctrlcomb;
        combined.val_RDratio=adjustedvalcomb;
        
        combined.minus=subrun1;
        combined.minus.combidx=ia;
        
        combined.plus=subrun2;
        combined.plus.combidx=ib;
        
        fc1=10.^(combined.RDratio(:,1)-3);
        fc2=10.^(combined.RDratio(:,2)-3);
        n1=subrun1.totalcount(ia(hasapt1));
        n2=subrun2.totalcount(ib(hasapt1));
        
        foldchange=fc2./fc1;
        p=(fc1.*n1+fc2.*n2)./(n1+n2);
        z=(fc1-fc2)./sqrt(p.*(1-p).*(1./n1+1./n2));

               
        isdiff=(abs(z)>-norminv(alpha/2)) & foldchange>1.5;
        combined.isdiff=isdiff;
        combined.foldchange=foldchange;
%         plot(combined.RDratio(combined.isdiff,1),combined.RDratio(combined.isdiff,2),'.','color',[0.4 0.4 0.9],'markersize',8)

        H=[];
        P=[];
        for h=1:length(RD1(:,1))
            [hyp,pva]=ttest2(RD1(h,:),RD2(h,:),'alpha',alpha);
            H(end+1)=hyp;
            P(end+1)=pva;
        end
%         [H,p]=ttest2(RD1,RD2,'alpha',alpha,'dim',2);
        
        combined.H=find(H(hasapt1));
        combined.P=P(hasapt1);
        
%         plot(combined.RDratio(combined.H,1),combined.RDratio(combined.H,2),'.','color',[0.9 0.4 0.7],'markersize',8)
        
        plot(adjustedctrlcomb([1 2 3 4 5],1),adjustedctrlcomb([1 2 3 4 5],2),'+','linewidth',2,'markersize',14,'color',[0.9 0.2 0.2])
%         plot(adjustedctrlcomb([ 2  4 5],1),adjustedctrlcomb([ 2  4 5],2),'+','linewidth',2,'markersize',14)
        plot(adjustedvalcomb(:,1),adjustedvalcomb(:,2),'o','linewidth',2,'markersize',14,'color',[0.9 0.7 0.2])
        
%         plot(x,fitcoeffs(1)+fitcoeffs(2)*x,'k--','linewidth',2)
        x=linspace(min(min(ctrlcomb([2 4 5],:))),max(max(ctrlcomb([2 4 5],:))));
        plot(x,x,'k--','linewidth',2)
        set(gca,'fontsize',22)
        set(gca,'linewidth',2)
    
        
        fc1=10.^(combined.ctrls_RDratio(:,1)-3);
        fc2=10.^(combined.ctrls_RDratio(:,2)-3);
        n1=subrun1.ctrlcounts;
        n2=subrun2.ctrlcounts;
%         
%         p=(fc1.*n1+fc2.*n2)./(n1+n2);
%         z=(fc1-fc2)./sqrt(p.*(1-p).*(1./n1+1./n2));
%                
%         ctrlsisdiff=abs(z)>-norminv(alpha/2);
%         combined.ctrlsisdiff=ctrlsisdiff;

        combined.totalcount=subrun1.totalcount(ia(hasapt1),:);
        combined.ctrlcounts=subrun1.ctrlcounts;
        combined.ctrlseqs=subrun1.ctrls.seqs;
        combined.valcounts=subrun1.valcounts;
        combined.valseqs=subrun1.val.seqs;

% xan_all_RSVIII.totalcount=xan_minus_RSVIII.totalcount;
% xan_all_RSVIII.ctrlcounts=xan_minus_RSVIII.ctrlcounts;
% xan_all_RSVIII.ctrlseqs=RSVIIItrim.ctrls.seqs;
% xan_all_RSVIII.valcounts=xan_minus_RSVIII.valcounts;
% xan_all_RSVIII.valseqs=RSVIIItrim.seqs(RSVIIItrim.val.ids);
        
    end
    
    function plotreplicates(~,minus,plus)
        RDrep1=[minus.adjustedRDratio(:,1);plus.adjustedRDratio(:,1)];
        RDrep2=[minus.adjustedRDratio(:,2);plus.adjustedRDratio(:,2)];

        crep1=[minus.adjustedRDratio_ctrls(:,1);plus.adjustedRDratio_ctrls(:,1)];
        crep2=[minus.adjustedRDratio_ctrls(:,2);plus.adjustedRDratio_ctrls(:,2)];

        vrep1=[minus.adjustedRDratio_val(:,1);plus.adjustedRDratio_val(:,1)];
        vrep2=[minus.adjustedRDratio_val(:,2);plus.adjustedRDratio_val(:,2)];

        scatter(RDrep1,RDrep2,70,'filled','markerfacealpha',0.2,'markerfacecolor',[0.1 0.1 0.1])
        hold on
        plot(crep1,crep2,'+','linewidth',2,'markersize',14)
        
        plot(vrep1,vrep2,'o','linewidth',2,'markersize',14) % comment if cdGI
        
        mdl=fitlm(RDrep1,RDrep2);
        rsq=mdl.Rsquared.Adjusted;
        fitcoeffs=[mdl.Coefficients.Estimate(1) mdl.Coefficients.Estimate(2)];
        x=linspace(min(RDrep1),max(RDrep2));

        plot(x,fitcoeffs(1)+fitcoeffs(2)*x,'k--','linewidth',2)


        text(1,2,sprintf('R^2=%0.2f',rsq),'fontsize',20)
        set(gca,'fontsize',24)
        set(gca,'linewidth',2)
        axis([.5 2.2 .5 2.2])
        xlabel('log_1_0(RNA/DNA) replicate 1')
        ylabel('log_1_0(RNA/DNA) replicate 2')
        box on
        legend('Library sequence','Control ribozymes','Validation candidates','1:1','location','southeast')
%         legend('Library sequence','Control ribozymes','1:1','location','southeast') % uncomment if cdGI
    end
    
    function [dimx,dimy]=findsubplotdim(~,ind)
        dimy=ceil(sqrt(ind));
        dimx=dimy;
        if (dimx*dimy)~=ind
            while dimy*(dimx)>ind
                if (dimx-1)*dimy>=ind
                    dimx=dimx-1;
                else
                    break;
                end
            end
        end
    end
    function [mdl, rsq]=plotlinfit(~,x,y)
    xval=linspace(min(x),max(x));
    mdl=fitlm(x,y);
    rsq=mdl.Rsquared.Adjusted;

    fitcoeffs=[mdl.Coefficients.Estimate(1) mdl.Coefficients.Estimate(2)];
    plot(xval,fitcoeffs(1)+fitcoeffs(2)*xval,'k--','linewidth',2)

    end
    
    function out=findfold(~,combined)
        out=combined;
        out.fold=10.^(combined.RDratio(:,2)-combined.RDratio(:,1));
        out.tapprox=1./sqrt(combined.totalcount);
        
    end
    function findMutants(data,parentseq,familycutoff,familynames)
        for k=1:length(parentseq)
            data.mutantfamily(k).dist=[];
            data.mutantfamily(k).seqidx=[];
            data.mutantfamily(k).seqs={};
            data.mutantfamily(k).counts=[];
            data.mutantfamily(k).name=familynames{k};
            data.mutantfamily(k).parentseq=parentseq{k};
            if isempty(familycutoff)
                familycutoff=5;
            end
            for j=1:length(data.seqs)
                d=seqpdist({data.seqs{j},parentseq{k}})*length(parentseq{k});
                if d<=familycutoff
                    data.mutantfamily(k).seqidx(end+1)=j;
                    data.mutantfamily(k).seqs{end+1}=data.seqs{j};
                    data.mutantfamily(k).dist(end+1)=d;
                    data.mutantfamily(k).counts(end+1)=data.counts(j,k);
                end
            end
        end
        
    end
    function [pairaligned,missense,insertion,deletion]=findMutScan(data,seqs,counts,library)
        if nargin<2
            library=1:length(data.mutantfamily);
        end        

        pairaligned=struct;
        
        for k=library
        p=data.mutantfamily(k).parentseq;
        % align with parent
        if nargin<2
            seqs=data.mutantfamily(k).seqs;
            counts=data.mutantfamily(k).counts;
        end
        insertion=zeros(length(p),4);
        deletion=zeros(length(p),1);
        missense=zeros(length(p),4);
        % s1=data.mutantfamily(k).seqs(sum(data.mutantfamily(k).compareTRPinter.countsminus,2)==max(sum(data.mutantfamily(k).compareTRPinter.countsminus,2)));
        DNA={'A','C','G','T','-'};
        pairaligned(k).seqs={};
        for i=1:length(seqs)
        % for i=1:1000
            s1=seqs(i);
            [s,a]=nwalign(p,s1,'Alphabet','NT');

            alignedseq=a{1}(3,:);
            pairaligned(k).seqs{end+1}=alignedseq;
            mismatch=find(a{1}(2,:)~='|');

            mutated2base=alignedseq(mismatch);
            parentbase=a{1}(1,mismatch);

            ins=find(parentbase=='-');
            del=find(mutated2base=='-');
            mis=find(parentbase~='-'&mutated2base~='-');

            if ~isempty(ins)
                if length(ins)<6

                for j=1:length(DNA)
                    base=find(mutated2base(ins)==DNA{j});
                    if ~isempty(base)
                        insertion(mismatch(ins),j)=insertion(mismatch(ins),j)+sum(counts(i));
                    end
                end
                end
            end



            if ~isempty(del)
                deletion(mismatch(del))=deletion(mismatch(del))+sum(counts(i));
            end

            if ~isempty(mis)
                if ~isempty(ins)

                for j=1:length(DNA)
                    base=find(mutated2base(mis)==DNA{j});

                    if ~isempty(base)
                            for q=1:length(base)
                                subtractlen=sum(ins<(base(q)));
                                missense(mismatch(mis)-subtractlen,j)=missense(mismatch(mis)-subtractlen,j)+sum(counts(i));
                            end

                    end
                end
                else
                    for j=1:length(DNA)
                        base=find(mutated2base(mis)==DNA{j});

                        if parentbase(mis)~=p(mismatch(mis))
                            fprintf('%s\n%s\n',p,a{1}(1,:))
                        end

                        if mutated2base(mis(base))==p(mismatch(mis(base)))
                            fprintf('%s\n%s\n%s\n-\n',p,a{1}(1,:),alignedseq)
                        end

                        if ~isempty(base)
                            missense(mismatch(mis(base)),j)=missense(mismatch(mis(base)),j)+sum(counts(i));

                            if j==2 & ismember(79,mismatch(mis(base)))
                            fprintf('%s\n%s\n%s\n-\n',p,a{1}(1,:),alignedseq)
                            end
                        end
                    end
                end
            end
        end
        

        missense=missense/sum(sum(counts));
        deletion=deletion/sum(sum(counts));
        insertion=insertion/sum(sum(counts));
%         
        
        setfig(strcat(data.mutantfamily(k).name,' mut scan'));clf
        hold on
        for q=1:length(missense(1,:))
            s=scatter(1:length(p),missense(:,q),'filled');
            s.MarkerFaceColor=s.CData;
        end
        s=scatter(1:length(p),deletion,'filled');
        s.MarkerFaceColor=s.CData;
        s=scatter(1:length(p),sum(insertion,2),'filled');
        s.MarkerFaceColor=s.CData;
        set(gca,'YScale','log')
%         ylim([10^-4.5,10^-2.5])
        
        xlim([0,length(p)+1])
        set(gca,'XTick',[0 1:length(p) length(p)+1])
        set(gca,'XTickLabel',{'',p(:),''})
        set(gca,'fontsize',14)
        set(gca,'linewidth',1.5)
        ylabel('mutation frequency')
        title(data.mutantfamily(k).name)
        legend({DNA{:},'+'})
        grid on
        end
    end
    
    
    function [pairaligned,missense,insertion,deletion]= findMutScanNoIndel(data,seqs,counts,library)
        if nargin<2
            
            library=1:length(data.mutantfamily);
        end
        pairaligned=struct;
        
        for k=library
            if nargin<2
                seqs=data.mutantfamily(k).seqs;
                counts=data.mutantfamily(k).counts;
            end
        p=data.mutantfamily(k).parentseq;
        % align with parent
        insertion=zeros(length(p),4);
        deletion=zeros(length(p),1);
        missense=zeros(length(p),4);
        % s1=data.mutantfamily(k).seqs(sum(data.mutantfamily(k).compareTRPinter.countsminus,2)==max(sum(data.mutantfamily(k).compareTRPinter.countsminus,2)));
        DNA={'A','C','G','T'};
        pairaligned(k).misseqs={};
        pairaligned(k).delseqs={};
        pairaligned(k).inseqs={};
        pairaligned(k).mispos={};
        for i=1:length(seqs)
        % for i=1:1000
            s1=seqs(i);
            [s,a]=nwalign(p,s1,'Alphabet','NT');

            alignedseq=a{1}(3,:);
            
            if length(s1{1})==length(p) & isempty(regexp(alignedseq,'-'));
                mispos=find(a{1}(2,:)~='|');
                pairaligned(k).misseqs{end+1}=alignedseq;
                pairaligned(k).mispos{end+1}=mispos;
            elseif length(s1{1})<length(p)
                pairaligned(k).delseqs{end+1}=alignedseq;
            elseif length(s1{1})>length(p);
                pairaligned(k).inseqs{end+1}=alignedseq;
            end
            
        end
        mismatsum=zeros(length(p),4);
        for j=1:length(DNA)
            mismat=zeros(length(pairaligned(k).misseqs),length(p));
            for i=1:length(pairaligned(k).misseqs)
                mut2base=pairaligned(k).misseqs{i}(pairaligned(k).mispos{i})==DNA{j};
                mismat(i,pairaligned(k).mispos{i}(mut2base))=sum(counts(i,:));
                mismat(i,pairaligned(k).mispos{i}(mut2base))=1;
            
            end
            mismatsum(:,j)=(sum(mismat)');
        end
            
        missense=mismatsum/sum(sum(counts));
        end
        
% 
%         missense=missense/sum(sum(counts));
%         deletion=deletion/sum(sum(counts));
%         insertion=insertion/sum(sum(counts));
% %         
%         
        setfig(strcat(data.mutantfamily(k).name,' missense scan'));clf
        hold on
        for q=1:length(missense(1,:))
            s=scatter(1:length(p),missense(:,q),'filled');
            s.MarkerFaceColor=s.CData;
        end
%         s=scatter(1:length(p),deletion,'filled');
%         s.MarkerFaceColor=s.CData;
%         s=scatter(1:length(p),sum(insertion,2),'filled');
%         s.MarkerFaceColor=s.CData;
        set(gca,'YScale','log')
%         ylim([10^-4.5,10^-2.5])
        
        xlim([0,length(p)+1])
        set(gca,'XTick',[0 1:length(p) length(p)+1])
        set(gca,'XTickLabel',{'',p(:),''})
        set(gca,'fontsize',14)
        set(gca,'linewidth',1.5)
        ylabel('mutation frequency')
        title(data.mutantfamily(k).name)
        legend({DNA{:},'+'})
        grid on
    end
    
    
    
    
   end
end
  