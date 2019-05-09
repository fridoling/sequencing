
# TO DO list for sequencing analysis
- general strategy to filter mutations found by isomut:
merge all the mutations, without any additional filtering after isomut
keep only those (chromosome / position) above a given threshold in at least one strain
keep those mutations also when they're below the threshold, to have a track of the mutations shared among different strains

which threshold shall we use?
we can try with: reads >= 3 (stricter than before), mut_freq >= .21, cleanliness >= .93, score >= .21
Paolo

- with this strategy, we could compare the frequency of mutations at 18C and 30C: those that decreased their frequency at 30C are likely to disturb cells when microtubules are correctly polymerized.
Paolo

- question: chromosomal duplication occurs before or after 'interesting' mutations?
if chrVIII duplication occurred before a mutation in that chr, we can find it at 100% frequency. otherwise, it can't be more than 50%.
at which frequency do we find chrVIII mutations in strains where chrVIII is duplicated?
the same question for other chromosomes.
Paolo
(Fridolin, can you produce a txt file with the median coverage per chromosome in every strain, to be included automatically in the R pipeline?)

- together with the coverage we should include in the R pipeline the data from the FACS analysis: who is diploid, who is haploid, who has heterogeneous ploidy.
Paolo

- are there chromosomes with different mutation frequency than the average?
Paolo

- some strains seem to have no mutations but intergenic or synonymous (C7 at 18C, for example).
who are them?
for this reason, we must also include intergenic/synonymous mutations in our analysis.
at least, we should include intergenic mutations and analyze more carefully the genes that surround them.
Paolo

- make cluster great again:
 move all the data, install all the software, write a pipeline to create bams 
DONE: it works and the reads in the new bams are the same as the old ones (according to picard’s CompareSAM)
run isomut on the cluster
Paolo and Fridolin

- the working place is shares.ifom: a new folder will be created into
 /ANC/evolution_tub2/tub2_401_TRIAL2_2018_mattia/11_NGS_analysis/
with everything
Paolo / Fridolin

## other things:
- we need a finalized version of the bam file preparation pipeline.
in the pipeline we're using now, we don't re-align around indels, we don’t have a saved copy of the output of each program (useful to find problems/errors/warnings) and we miss a final check of bam file consistency.
once implemented these steps, we should run it (on the cluster) on all the raw data we already have, and repeat all the isomut analysis on the new bam files, to have them consistent.
Fridolin

- we never found chromosomal rearrangements: are we missing them? which software shall we use to spot them?
Fridolin

- include a gene enrichment / gene ontology analysis of the mutations
Fridolin

- compute the expected number of mutations (from Petrov, PNAS, and other papers)
Paolo/Andrea/everybody

- have a finalized list of mutations present in the ancestors, to look for them in the evolved strains (using samtools, not isomut)
Fridolin?

- What reference should we use? I looked at the methods part of the Petrov paper (Zhu et al. 2014, PNAS). They create their own reference genome based on shared mutations with respect to S288C and then repeat the alignment of all their samples using that modified reference “to fully eliminate confounding influences from these fixed differences that may affect mapping accuracy”. Something that we might also think about. Aside from that  there is the paper presenting the full W303 sequence (Matheson et al. 2017, Genes, Genomes, Genetics) that also lists all the variants with respect to S288C. The sequence itself can be found here: https://www.ncbi.nlm.nih.gov/Traces/wgs/LYZE01

- about intergenic mutations: for example, in C7 at 18C there are three, all at high frequency.
 we could take a clone out the well, check that the three mutations are there, cross the clone with a wild type, and analyze the cold-sensitivity of the segregants.
we could leave this step after the next sequencing and analysis, once we trust them more.
Elena, after the next NGS

