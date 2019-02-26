
Rscript plot2DO_setup.R

Rscript plot2DO.R --file=sample_BAM_files/CC_yeast_GSM2561059.5M_reads.bed \
                  --reference=Plus1 --minLength=20 --maxLength=200

Rscript plot2DO.R --file=sample_BAM_files/yeast_SRR2045610.5M_reads.bam

Rscript plot2DO.R --file=sample_BAM_files/yeast_SRR2045610.5M_reads.bam \
                  --type=dyads --reference=Plus1 --colorScaleMax=0.15

Rscript plot2DO.R --file=sample_BAM_files/yeast_SRR2045610.5M_reads.bam \
                  --sites=annotations/Yeast_tRNA_genes.bed \
                  --align=fivePrime --siteLabel=tRNA --simplifyPlot=on

Rscript plot2DO.R --file=sample_BAM_files/yeast_SRR2045610.5M_reads.bam \
                  --sites=annotations/Yeast_ARS_locations.bed \
                  --align=center --siteLabel=ARS --simplifyPlot=on

Rscript plot2DO.R --file=sample_BAM_files/yeast_SRR2045610.5M_reads.bam \
                  --sites=annotations/Yeast_centromeres.bed \
                  --align=threePrime --siteLabel=centromere --simplifyPlot=on

Rscript plot2DO.R --file=sample_BAM_files/fly_SRR2038265.5M_reads.bam \
                  --genome=dm3 --reference=TSS --simplifyPlot=on

Rscript plot2DO.R --file=sample_BAM_files/worm_SRR3289717.5M_reads.bam \
                  --genome=ce10 --reference=TSS --simplifyPlot=on

Rscript plot2DO.R --file=sample_BAM_files/mouse_SRR572708.5M_reads.bam \
                  --genome=mm10 --reference=TSS --simplifyPlot=on

Rscript plot2DO.R --file=sample_BAM_files/human_SRR1781839.5M_reads.bam \
                  --genome=hg19 --reference=TSS --simplifyPlot=on

Rscript plot2DO.R --file=sample_BAM_files/yeast_50U_MNase_SRR3649301.5M_reads.bam \
                  --sites=annotations/Yeast_tRNA_genes.bed \
                  --siteLabel=tRNA --colorScaleMax=0.03 \
                  --upstream=100 --downstream=100 --squeezePlot=on

Rscript plot2DO.R --file=sample_BAM_files/yeast_100U_MNase_SRR3649296.5M_reads.bam \
                  --sites=annotations/Yeast_tRNA_genes.bed \
                  --siteLabel=tRNA --colorScaleMax=0.03 \
                  --upstream=100 --downstream=100 --squeezePlot=on

Rscript plot2DO.R --file=sample_BAM_files/yeast_200U_MNase_SRR3649297.5M_reads.bam \
                  --sites=annotations/Yeast_tRNA_genes.bed \
                  --siteLabel=tRNA --colorScaleMax=0.03 \
                  --upstream=100 --downstream=100 --squeezePlot=on

Rscript plot2DO.R --file=sample_BAM_files/yeast_300U_MNase_SRR3649298.5M_reads.bam \
                  --sites=annotations/Yeast_tRNA_genes.bed \
                  --siteLabel=tRNA --colorScaleMax=0.03 \
                  --upstream=100 --downstream=100 --squeezePlot=on

Rscript plot2DO.R --file=sample_BAM_files/yeast_400U_MNase_SRR3649299.5M_reads.bam \
                  --sites=annotations/Yeast_tRNA_genes.bed \
                  --siteLabel=tRNA --colorScaleMax=0.03 \
                  --upstream=100 --downstream=100 --squeezePlot=on

Rscript plot2DO.R --file=sample_BAM_files/yeast_50U_MNase_SRR3649301.5M_reads.bam \
                  --reference=Plus1 --colorScaleMax=0.03 \
                  --upstream=100 --downstream=100 --squeezePlot=on

Rscript plot2DO.R --file=sample_BAM_files/yeast_100U_MNase_SRR3649296.5M_reads.bam \
                  --reference=Plus1 --colorScaleMax=0.03 \
                  --upstream=100 --downstream=100 --squeezePlot=on

Rscript plot2DO.R --file=sample_BAM_files/yeast_200U_MNase_SRR3649297.5M_reads.bam \
                  --reference=Plus1 --colorScaleMax=0.03 \
                  --upstream=100 --downstream=100 --squeezePlot=on

Rscript plot2DO.R --file=sample_BAM_files/yeast_300U_MNase_SRR3649298.5M_reads.bam \
                  --reference=Plus1 --colorScaleMax=0.03 \
                  --upstream=100 --downstream=100 --squeezePlot=on

Rscript plot2DO.R --file=sample_BAM_files/yeast_400U_MNase_SRR3649299.5M_reads.bam \
                  --reference=Plus1 --colorScaleMax=0.03 \
                  --upstream=100 --downstream=100 --squeezePlot=on
