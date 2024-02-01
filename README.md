# A collection of scripting software used in the stickleback project

Marine Assembly: From Pacbio, Nanopore, 10x Genomics, to HiC
Lastz workflow: axt to chain, to netted chains

## GATK SNPs/INDELs calling work flow

#### Reference

<https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels>

### Step 0: Environment setup, (see dockerfile)

<https://github.com/edotau/sticklebackCipher/blob/db6adda15296449debd2b224d2ac975da4cec9c1/gatk/Dockerfile#L1-L18>
TODO: conda environment might be the better option to set up env and third party software

### Step 1: Trim Fastqs and align to reference

<https://github.com/edotau/sticklebackCipher/blob/8906a84f3471bbd2c7aa741e3c92dd2d37372616/gatk/src/alignReads.bash#L43-L66>

### Step 2: Process BAMs and call raw genotype variance

<https://github.com/edotau/sticklebackCipher/blob/8906a84f3471bbd2c7aa741e3c92dd2d37372616/gatk/src/gatkfmtBam.bash#L37-L49>

### Step 3: Combine raw genotypes into one population vcf

<https://github.com/edotau/sticklebackCipher/blob/8906a84f3471bbd2c7aa741e3c92dd2d37372616/gatk/src/mergeGvcf.bash#L27-L41>

### Step 4: Initial filtering of raw snp and indel calls

<https://github.com/edotau/sticklebackCipher/blob/acc009d83ed472bb6f661e2bf0685c73c2d9fddb/gatk/src/snpIndelDiscovery.bash#L18-L27>

### Step 5: Recalibrate initial alignments and recall variance using genotype vcf

<https://github.com/edotau/sticklebackCipher/blob/acc009d83ed472bb6f661e2bf0685c73c2d9fddb/gatk/src/baseRecalibratorBQSR.bash#L29-L50>

### Step 6: Combine population genotypes into one vcf

<https://github.com/edotau/sticklebackCipher/blob/acc009d83ed472bb6f661e2bf0685c73c2d9fddb/gatk/src/mergeGvcf.bash#L27-L41>

### Step 7: The variant recalibration step fits a Gaussian mixture model to the contextual annotations given to each variant

<https://github.com/edotau/sticklebackCipher/blob/acc009d83ed472bb6f661e2bf0685c73c2d9fddb/gatk/src/gatkGenotypeVariants.bash#L35-L62>
