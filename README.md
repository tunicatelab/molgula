# _Molgula occulta_, _M. oculata_ and hybrid analysis (Lowe et al. 2020, Fodor et al. 2021)
_M. occulta_ and _M. oculata_ are tunicate species found on the Northwestern coast of France, in Roscoff. They look very similar in their adult stages, however have alternate body plans in their larval stage. _M. oculata_ goes through typical tunicate development, with a tadpole larvae. _M. occulta_, on the hand forgoes developing a tail and goes into metamorphosis shortly after hatching. These species are able to be cross-fertilized producing a hybrid with half a tail. Using RNA-seq, we sequenced and examined three developmental time points for _M. occulta, M. oculata_ and their hybrid (Lowe et al. 2014). Their genomes were sequenced and assembled in parallel to a third species, _M. occidentalis_ (Stolfi et al. 2014).

Genomes from [Aniseed](https://aniseed.cnrs.fr) were scaffolded using [redundans](https://github.com/lpryszcz/redundans) with reads from the original assemblies, and _M. oculata_ being used as reference for _M. occulta_ because _M. oculata_ assembly is less fragmented and closely related. 

Transcriptomes were assembled using [Trinity genome guided command](https://github.com/trinityrnaseq/RagonInst_Sept2017_Workshop/wiki/genome_guided_trinity) this was done to improve the gene models that are currently found on ANISEED. Next, transcripts were converted into [supertranscripts](https://github.com/trinityrnaseq/trinityrnaseq/wiki/SuperTranscripts) to give a gene view to aid in cross-species comparisons. Transcripts were then filtered for short (less than 100 amino acids) fragments using [TransDecoder](https://github.com/TransDecoder/TransDecoder). Finally, transcripts were scaffolded with TransPS](https://bioinformatics.cs.vt.edu/zhanglab/software/transps/) using [*C. robusta*](https://www.aniseed.cnrs.fr/aniseed/download/download_data) protein sequences as a reference. 
Genome and transcriptome assemblies can be found on [open science framework](https://osf.io/mj3r7/)

_M. occulta_ and _M. oculata_ reads were mapped to their respective transcriptomes using [bowtie2](https://github.com/BenLangmead/bowtie2) and [eXpress](https://pachterlab.github.io/eXpress/overview.html) streaming pipeline. Reads from _occulta X oculata_ hybrids were mapped onto both parent transcriptomes in the same manner in order to quantify allele-specific expression. 

[EdgeR](https://www.bioconductor.org/packages/devel/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf) was used to identify differentially expressed genes. 

REFERENCES:

Lowe EK, Swalla BJ, Brown CT. Evaluating a lightweight transcriptome assembly pipeline on two closely related ascidian species. PeerJ PrePrints; 2014 Sep 18.

Stolfi A, Lowe EK, Racioppi C, Ristoratore F, Brown CT, Swalla BJ, Christiaen L. Divergent mechanisms regulate conserved cardiopharyngeal development and gene expression in distantly related ascidians. elife. 2014 Sep 10;3:e03728.

Lowe EK, Racioppi C, Peyriéras N, Ristoratore F, Christiaen L, Swalla BJ, Stolfi A. A cis‐regulatory change underlying the motor neuron‐specific loss of Ebf expression in immotile tunicate larvae. Evolution & Development. 2020 Dec 23:e12364.
