"""
Genome, genetic map and demographic model definitions for humans.
"""
import math

import msprime

import stdpopsim.models as models
import stdpopsim.genomes as genomes
import stdpopsim.genetic_maps as genetic_maps


###########################################################
#
# Genetic maps
#
###########################################################


class HapmapII_GRCh37(genetic_maps.GeneticMap):
    """
    The Phase II HapMap Genetic map (lifted over to GRCh37) used in
    1000 Genomes. Please see the `README
    <ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/working/20110106_recombination_hotspots/README_hapmapII_GRCh37_map>`_
    for more details.
    """
    url = (
        "http://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/working/"
        "20110106_recombination_hotspots/"
        "HapmapII_GRCh37_RecombinationHotspots.tar.gz")
    file_pattern = "genetic_map_GRCh37_{name}.txt"


genetic_maps.register_genetic_map(HapmapII_GRCh37())


###########################################################
#
# Genome definition
#
###########################################################

# List of chromosomes.

# FIXME: add mean mutation rate data to this table.
# Name  Length  mean_recombination_rate mean_mutation_rate

# length information can be found here
# <http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/chromInfo.txt.gz>

# mean_recombination_rate was computed across all windows of the GRCh37 genetic map
# <ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/working/20110106_recombination_hotspots>
_chromosome_data = """\
chr1 	 249250621 	 1.1485597641285933e-08
chr2 	 243199373 	 1.1054289277533446e-08
chr3 	 198022430 	 1.1279585624662551e-08
chr4 	 191154276 	 1.1231162636001008e-08
chr5 	 180915260 	 1.1280936570022824e-08
chr6 	 171115067 	 1.1222852661225285e-08
chr7 	 159138663 	 1.1764614397655721e-08
chr8 	 146364022 	 1.1478465778920576e-08
chr9 	 141213431 	 1.1780701596308656e-08
chr10 	 135534747 	 1.3365134257075317e-08
chr11 	 135006516 	 1.1719334320833283e-08
chr12 	 133851895 	 1.305017186986983e-08
chr13 	 115169878 	 1.0914860554958317e-08
chr14 	 107349540 	 1.119730771394731e-08
chr15 	 102531392 	 1.3835785893339787e-08
chr16 	 90354753 	 1.4834607113882717e-08
chr17 	 81195210 	 1.582489036239487e-08
chr18 	 78077248 	 1.5075956950023575e-08
chr19 	 59128983 	 1.8220141872466202e-08
chr20 	 63025520 	 1.7178269031631664e-08
chr21 	 48129895 	 1.3045214034879191e-08
chr22 	 51304566 	 1.4445022767788226e-08
chrX 	 155270560 	 1.164662223273842e-08
chrY 	 59373566 	 0.0
"""

_chromosomes = []
for line in _chromosome_data.splitlines():
    name, length, mean_rr = line.split()[:3]
    _chromosomes.append(genomes.Chromosome(
        name=name, length=int(length),
        default_mutation_rate=1e-8,  # WRONG!,
        default_recombination_rate=float(mean_rr)))


#: :class:`stdpopsim.Genome` definition for humans.
genome = genomes.Genome(
    species="homo_sapiens",
    chromosomes=_chromosomes,
    default_genetic_map=HapmapII_GRCh37.name)


###########################################################
#
# Demographic models
#
###########################################################


class GutenkunstThreePopOutOfAfrica(models.Model):
    """
    The three population Out-of-Africa model from Gutenkunst et al.

    .. todo:: document this model, including the original publications
        and clear information about what the different population indexes
        mean.

    """

    def __init__(self):
        super().__init__()
        # First we set out the maximum likelihood values of the various parameters
        # given in Table 1.
        N_A = 7300
        N_B = 2100
        N_AF = 12300
        N_EU0 = 1000
        N_AS0 = 510
        # Times are provided in years, so we convert into generations.
        generation_time = 25
        T_AF = 220e3 / generation_time
        T_B = 140e3 / generation_time
        T_EU_AS = 21.2e3 / generation_time
        # We need to work out the starting (diploid) population sizes based on
        # the growth rates provided for these two populations
        r_EU = 0.004
        r_AS = 0.0055
        N_EU = N_EU0 / math.exp(-r_EU * T_EU_AS)
        N_AS = N_AS0 / math.exp(-r_AS * T_EU_AS)
        # Migration rates during the various epochs.
        m_AF_B = 25e-5
        m_AF_EU = 3e-5
        m_AF_AS = 1.9e-5
        m_EU_AS = 9.6e-5
        # Population IDs correspond to their indexes in the population
        # configuration array. Therefore, we have 0=YRI, 1=CEU and 2=CHB
        # initially.
        self.population_configurations = [
            msprime.PopulationConfiguration(initial_size=N_AF),
            msprime.PopulationConfiguration(initial_size=N_EU, growth_rate=r_EU),
            msprime.PopulationConfiguration(initial_size=N_AS, growth_rate=r_AS)
        ]
        self.migration_matrix = [
            [      0, m_AF_EU, m_AF_AS],  # noqa
            [m_AF_EU,       0, m_EU_AS],  # noqa
            [m_AF_AS, m_EU_AS,       0],  # noqa
        ]
        self.demographic_events = [
            # CEU and CHB merge into B with rate changes at T_EU_AS
            msprime.MassMigration(
                time=T_EU_AS, source=2, destination=1, proportion=1.0),
            msprime.MigrationRateChange(time=T_EU_AS, rate=0),
            msprime.MigrationRateChange(
                time=T_EU_AS, rate=m_AF_B, matrix_index=(0, 1)),
            msprime.MigrationRateChange(
                time=T_EU_AS, rate=m_AF_B, matrix_index=(1, 0)),
            msprime.PopulationParametersChange(
                time=T_EU_AS, initial_size=N_B, growth_rate=0, population_id=1),
            # Population B merges into YRI at T_B
            msprime.MassMigration(
                time=T_B, source=1, destination=0, proportion=1.0),
            # Size changes to N_A at T_AF
            msprime.PopulationParametersChange(
                time=T_AF, initial_size=N_A, population_id=0)
        ]


class TennessenEuropean(models.Model):
    """
    The model is derived from the Tennesen et al.
    `analysis <https://doi.org/10.1126/science.1219240>`_  of the jSFS from
    European Americans and African Americans.

    .. todo:: document this model, including the original publications
        and clear information about what the different population indexes
        mean.

    """
    def __init__(self):
        super().__init__()
        # Population sizes
        N_A = 7310
        N_AF = 14474
        N_B = 1861
        N_EU1 = 9475
        # Times
        generation_time = 25
        T_AF = 148000 / generation_time
        T_B = 51000 / generation_time
        T_EU0 = 23000 / generation_time
        T_EU1 = 5115 / generation_time
        # Rates and Present Ne
        r_EU0 = 0.00307
        r_EU1 = 0.0195
        N_EU = N_EU1 / math.exp(-r_EU1 * T_EU1)
        self.population_configurations = [
            msprime.PopulationConfiguration(initial_size=N_EU, growth_rate=r_EU1)
        ]
        self.demographic_events = [
            msprime.PopulationParametersChange(
                time=T_EU1, initial_size=N_EU1, growth_rate=r_EU0, population_id=0),
            msprime.PopulationParametersChange(
                time=T_EU0, initial_size=N_B, growth_rate=0, population_id=0),
            msprime.PopulationParametersChange(
                time=T_B, initial_size=N_AF, growth_rate=0, population_id=0),
            msprime.PopulationParametersChange(
                time=T_AF, initial_size=N_A, growth_rate=0, population_id=0)
        ]


class GazaveTwoPopOutOfAfrica(models.Model):
    """
    docs
    
    """
    def __init__(self):
        super().__init__()
        generation_time = 25

        # 220kya:
        # African population constant with Ne~7300
        N_A = 7310
        
        # 148kya:
        # instantaneous growth to Ne~14000
        T_AF = 148e3 / generation_time
        N_AF = 14474
        
        N6_EU = 13143
        
        # 118kya:
        # non-AFR pops migrate OOA; bottlenecks to Ne~1800
        # migration between AFR occurs
        N_B = 1861
        T5 = 118e3 / generation_time
        T4 = T5
        m_AF_B = 15e-5
        N5_EU = 62
        N4_EU = N6_EU
        
        # 18kya:
        # 2nd EUR bottlenecks to Ne~1000 & starts growing with rate 0.307%
        # migration rate slows between AFR-EUR
        N_EU0 = 1032
        T3 = 18e3 / generation_time
        T2 = T3
    #     m_AF_EU = 2.5e-5
        r_EU0 = 0 # 0.00307
    #     N2_EU = 15829 # N_EU0 / math.exp(-r_EU0 * T_EU_B)
        N2_EU = 16178
        N2_AF = 26682
        N3_EU = 2020
        # 5.1kya:
        # explosive growth in both AFR & EUR
    
        # Chen 2015
    #     T1_EU = 7.26e3 / generation_time 
        T1_EU = 4.95e3 / generation_time
        T1_AF = 10.01e3 / generation_time
    #     r_EU = 0.0149
        r_EU = 0.022
        r_AF = 0.00735
    #     N1_EU = 1.2e6 # N_EU1 / math.exp(-r_EU * T_EG)
        N1_EU = 1.261e6
        m_EG = 0
        N1_AF = 5.062e5 # N_AF / math.exp(-r_AF * T_EG)
        
        # Population IDs correspond to their indexes in the population
        # configuration array. Therefore, we have 0=YRI, 1=CEU initially.
        self.population_configurations = [
            msprime.PopulationConfiguration(
                sample_size=n_afr, initial_size=N1_AF, growth_rate=r_AF),
            msprime.PopulationConfiguration(
                sample_size=n_eur, initial_size=N1_EU, growth_rate=r_EU)#,
        ]
    
        # up to 5.1kya, no migration
        self.migration_matrix = [
            [0, 0],
            [0, 0],
        ]
        
        self.demographic_events = [
            # at 5.1kya, change to slow growth rate in EUR & stop growth in AFR;
            # add migration rate
            msprime.PopulationParametersChange(
                time=T1_EU, growth_rate=0, initial_size=N2_EU, population_id=1),
            msprime.PopulationParametersChange(
                time=T1_AF, growth_rate=0, initial_size=N2_AF, population_id=0),
            
            # at 18kya, bottleneck + instantaneous recovery to smaller Ne
            msprime.PopulationParametersChange(
                time=T2, initial_size=N3_EU, population_id=1),
            msprime.PopulationParametersChange(
                time=T3, initial_size=N4_EU, population_id=1),
            # at 118kya, bottleneck + instantaneous recovery to same Ne
            msprime.PopulationParametersChange(
                time=T4, initial_size=N5_EU, population_id=1),
            msprime.PopulationParametersChange(
                time=T5, initial_size=N6_EU, population_id=1),
            # At 148kya, instantaneous growth in AFR
            msprime.PopulationParametersChange(
                time=T_AF, initial_size=N_A, population_id=0)
        ]
    


class ChenTwoPopOutOfAfrica(models.Model):
    """
    docs
    
    """
    def __init__(self):
        super().__init__()
        generation_time = 25

        # 220kya:
        # African population constant with Ne~7300
        N_A = 7310
        
        # 148kya:
        # instantaneous growth to Ne~14000
        T_AF = 148e3 / generation_time
        N_AF = 14474
        
        N6_EU = 13143
        
        # 118kya:
        # non-AFR pops migrate OOA; bottlenecks to Ne~1800
        # migration between AFR occurs
        N_B = 1861
        T5 = 118e3 / generation_time
        T4 = T5
        m_AF_B = 15e-5
        N5_EU = 62
        N4_EU = N6_EU
        
        # 18kya:
        # 2nd EUR bottlenecks to Ne~1000 & starts growing with rate 0.307%
        # migration rate slows between AFR-EUR
        N_EU0 = 1032
        T3 = 18e3 / generation_time
        T2 = T3
    #     m_AF_EU = 2.5e-5
        r_EU0 = 0 # 0.00307
    #     N2_EU = 15829 # N_EU0 / math.exp(-r_EU0 * T_EU_B)
        N2_EU = 16178
        N2_AF = 26682
        N3_EU = 2020
        # 5.1kya:
        # explosive growth in both AFR & EUR
    
        # Chen 2015
    #     T1_EU = 7.26e3 / generation_time 
        T1_EU = 4.95e3 / generation_time
        T1_AF = 10.01e3 / generation_time
    #     r_EU = 0.0149
        r_EU = 0.022
        r_AF = 0.00735
    #     N1_EU = 1.2e6 # N_EU1 / math.exp(-r_EU * T_EG)
        N1_EU = 1.261e6
        m_EG = 0
        N1_AF = 5.062e5 # N_AF / math.exp(-r_AF * T_EG)
        
        # Population IDs correspond to their indexes in the population
        # configuration array. Therefore, we have 0=YRI, 1=CEU initially.
        self.population_configurations = [
            msprime.PopulationConfiguration(
                sample_size=n_afr, initial_size=N1_AF, growth_rate=r_AF),
            msprime.PopulationConfiguration(
                sample_size=n_eur, initial_size=N1_EU, growth_rate=r_EU)#,
        ]
    
        # up to 5.1kya, no migration
        self.migration_matrix = [
            [0, 0],
            [0, 0],
        ]
        
        self.demographic_events = [
            # at 5.1kya, change to slow growth rate in EUR & stop growth in AFR;
            # add migration rate
            msprime.PopulationParametersChange(
                time=T1_EU, growth_rate=0, initial_size=N2_EU, population_id=1),
            msprime.PopulationParametersChange(
                time=T1_AF, growth_rate=0, initial_size=N2_AF, population_id=0),
            
            # at 18kya, bottleneck + instantaneous recovery to smaller Ne
            msprime.PopulationParametersChange(
                time=T2, initial_size=N3_EU, population_id=1),
            msprime.PopulationParametersChange(
                time=T3, initial_size=N4_EU, population_id=1),
            # at 118kya, bottleneck + instantaneous recovery to same Ne
            msprime.PopulationParametersChange(
                time=T4, initial_size=N5_EU, population_id=1),
            msprime.PopulationParametersChange(
                time=T5, initial_size=N6_EU, population_id=1),
            # At 148kya, instantaneous growth in AFR
            msprime.PopulationParametersChange(
                time=T_AF, initial_size=N_A, population_id=0)
        ]
    