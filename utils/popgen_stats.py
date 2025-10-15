"""
Pure Python population genetics statistics calculations
No external GenePop binary required
"""
import numpy as np
from scipy import stats
from collections import Counter, defaultdict


def parse_genotype(genotype_list):
    """
    BioPython already parses genotypes into list of tuples like [(1, 1), (2, 2)]
    Just return them as-is
    """
    if not genotype_list:
        return None
    return genotype_list


def calculate_allele_frequencies(record):
    """Calculate allele frequencies for each locus in each population"""
    results = {}

    for pop_idx, population in enumerate(record.populations):
        pop_name = f'Population_{pop_idx + 1}'
        results[pop_name] = {}

        for locus_idx, locus_name in enumerate(record.loci_list):
            allele_counts = Counter()
            total_alleles = 0

            # BioPython format: each individual is (name, [genotype_tuples])
            for individual_name, genotypes in population:
                if genotypes and locus_idx < len(genotypes) and genotypes[locus_idx]:
                    allele1, allele2 = genotypes[locus_idx]
                    if allele1 and allele1 != 0:  # 0 means missing
                        allele_counts[allele1] += 1
                        total_alleles += 1
                    if allele2 and allele2 != 0:
                        allele_counts[allele2] += 1
                        total_alleles += 1

            # Calculate frequencies
            frequencies = {}
            if total_alleles > 0:
                for allele, count in allele_counts.items():
                    frequencies[str(allele)] = round(count / total_alleles, 4)

            results[pop_name][locus_name] = {
                'frequencies': frequencies,
                'sample_size': total_alleles // 2  # diploid
            }

    return results


def calculate_heterozygosity(record):
    """Calculate observed and expected heterozygosity"""
    results = {}

    for pop_idx, population in enumerate(record.populations):
        pop_name = f'Population_{pop_idx + 1}'
        results[pop_name] = {}

        for locus_idx, locus_name in enumerate(record.loci_list):
            allele_counts = Counter()
            heterozygotes = 0
            total_individuals = 0

            for individual_name, genotypes in population:
                if genotypes and locus_idx < len(genotypes) and genotypes[locus_idx]:
                    allele1, allele2 = genotypes[locus_idx]
                    if allele1 and allele1 != 0 and allele2 and allele2 != 0:
                        allele_counts[str(allele1)] += 1
                        allele_counts[str(allele2)] += 1
                        total_individuals += 1

                        if allele1 != allele2:
                            heterozygotes += 1

            # Observed heterozygosity
            Ho = heterozygotes / total_individuals if total_individuals > 0 else 0

            # Expected heterozygosity (gene diversity)
            total_alleles = sum(allele_counts.values())
            He = 0
            if total_alleles > 0:
                frequencies = {allele: count / total_alleles for allele, count in allele_counts.items()}
                He = 1 - sum(freq ** 2 for freq in frequencies.values())

            results[pop_name][locus_name] = {
                'observed': round(Ho, 4),
                'expected': round(He, 4),
                'sample_size': total_individuals
            }

    return results


def hardy_weinberg_test(record):
    """Test Hardy-Weinberg equilibrium using chi-square test"""
    results = {}

    for pop_idx, population in enumerate(record.populations):
        pop_name = f'Population_{pop_idx + 1}'
        results[pop_name] = {}

        for locus_idx, locus_name in enumerate(record.loci_list):
            # Count genotypes
            genotype_counts = Counter()

            for individual_name, genotypes in population:
                if genotypes and locus_idx < len(genotypes) and genotypes[locus_idx]:
                    allele1, allele2 = genotypes[locus_idx]
                    if allele1 and allele1 != 0 and allele2 and allele2 != 0:
                        # Sort alleles to make genotype consistent
                        gt = tuple(sorted([str(allele1), str(allele2)]))
                        genotype_counts[gt] += 1

            # Calculate allele frequencies
            allele_counts = Counter()
            total = sum(genotype_counts.values())

            for (a1, a2), count in genotype_counts.items():
                allele_counts[a1] += count
                allele_counts[a2] += count

            total_alleles = sum(allele_counts.values())

            if total_alleles == 0 or len(allele_counts) < 2:
                results[pop_name][locus_name] = {'status': 'insufficient_data'}
                continue

            # Calculate expected genotype frequencies
            allele_freqs = {allele: count / total_alleles for allele, count in allele_counts.items()}

            # Calculate chi-square
            chi_square = 0
            expected_counts = {}

            for (a1, a2), observed in genotype_counts.items():
                if a1 == a2:  # Homozygote
                    expected = total * (allele_freqs[a1] ** 2)
                else:  # Heterozygote
                    expected = total * 2 * allele_freqs[a1] * allele_freqs[a2]

                expected_counts[f'{a1}/{a2}'] = round(expected, 2)

                if expected > 0:
                    chi_square += ((observed - expected) ** 2) / expected

            # Degrees of freedom = number of genotypes - number of alleles
            df = len(genotype_counts) - len(allele_counts)

            if df > 0:
                p_value = 1 - stats.chi2.cdf(chi_square, df)
                results[pop_name][locus_name] = {
                    'chi_square': round(chi_square, 4),
                    'p_value': round(p_value, 4),
                    'df': df,
                    'equilibrium': 'yes' if p_value > 0.05 else 'no',
                    'sample_size': total
                }
            else:
                results[pop_name][locus_name] = {'status': 'insufficient_variation'}

    return results


def calculate_fst(record):
    """Calculate Fst (fixation index) between populations"""
    if len(record.populations) < 2:
        return {'error': 'Need at least 2 populations'}

    results = {}

    for locus_idx, locus_name in enumerate(record.loci_list):
        # Calculate allele frequencies for each population
        pop_freqs = []
        pop_sizes = []

        for population in record.populations:
            allele_counts = Counter()
            total = 0

            for individual_name, genotypes in population:
                if genotypes and locus_idx < len(genotypes) and genotypes[locus_idx]:
                    allele1, allele2 = genotypes[locus_idx]
                    if allele1 and allele1 != 0:
                        allele_counts[str(allele1)] += 1
                        total += 1
                    if allele2 and allele2 != 0:
                        allele_counts[str(allele2)] += 1
                        total += 1

            if total > 0:
                pop_freqs.append({allele: count / total for allele, count in allele_counts.items()})
                pop_sizes.append(total // 2)
            else:
                pop_freqs.append({})
                pop_sizes.append(0)

        # Get all unique alleles
        all_alleles = set()
        for freqs in pop_freqs:
            all_alleles.update(freqs.keys())

        if len(all_alleles) == 0:
            continue

        # Calculate Fst using Weir & Cockerham method
        total_samples = sum(pop_sizes)

        Ht = 0  # Total heterozygosity
        Hs = 0  # Average within-population heterozygosity

        for allele in all_alleles:
            # Overall allele frequency
            p_mean = sum(pop_freqs[i].get(allele, 0) * pop_sizes[i] for i in range(len(pop_freqs))) / total_samples if total_samples > 0 else 0

            Ht += p_mean * (1 - p_mean)

            # Average within-population
            for i, freqs in enumerate(pop_freqs):
                p = freqs.get(allele, 0)
                Hs += (p * (1 - p) * pop_sizes[i]) / total_samples if total_samples > 0 else 0

        Fst = (Ht - Hs) / Ht if Ht > 0 else 0

        results[locus_name] = {
            'Fst': round(Fst, 4),
            'Ht': round(Ht, 4),
            'Hs': round(Hs, 4)
        }

    # Calculate average across loci
    if results:
        avg_fst = sum(r['Fst'] for r in results.values()) / len(results)
        results['average'] = round(avg_fst, 4)

    return results
