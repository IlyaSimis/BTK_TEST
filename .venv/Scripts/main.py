import sys
import random

def parse_vcf(vcf_file):
    """Парсит VCF-файл и возвращает словарь с вариантами."""

    variants = {}
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            cols = line.strip().split('\t')
            chrom = cols[0]
            pos = int(cols[1])
            ref = cols[3]
            alt = cols[4]
            samples = cols[len(cols):]
            if chrom not in variants:
                variants[chrom] = []
            variants[chrom].append((pos, ref, alt, samples))
    return variants


def parse_fasta(fasta_file):
    """Парсит FASTA-файл и возвращает словарь с референсным геномом."""

    sequences = {}
    current_chrom = None
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                current_chrom = line[1:].strip()
                sequences[current_chrom] = []
            else:
                sequences[current_chrom].append(line.strip())
    for chrom in sequences:
        sequences[chrom] = ''.join(sequences[chrom])
    return sequences


def generate_consensus(variants, reference, length, count, frequency_threshold):
    """Генерирует консенсусные последовательности."""

    consensus_sequences = []
    chrom_count = len(reference)
    per_chrom_count = count // chrom_count
    extra = count % chrom_count

    for i, (chrom, sequence) in enumerate(reference.items()):
        if chrom not in variants:
            continue

        chrom_variants = variants[chrom]
        chrom_length = len(sequence)
        chrom_positions = []
        current_count = per_chrom_count + (1 if i < extra else 0)

        while len(chrom_positions) < current_count:
            start = random.randint(0, chrom_length - length)
            if start in chrom_positions:
                continue
            chrom_positions.append(start)

            end = start + length
            consensus = list(sequence[start:end])
            for pos, ref, alt, samples in chrom_variants:
                if start <= pos < end:
                    sample_alleles = [s.split(':')[0] for s in samples if s]
                    if len(sample_alleles) > 0:  # Проверяем, что список не пуст
                        allele_freq = sum(1 for a in sample_alleles if a == '1') / len(sample_alleles)
                        if allele_freq >= frequency_threshold:
                            idx = pos - start
                            consensus[idx] = alt

            consensus_sequences.append(f'>{chrom}_{start}_{end}\n{"".join(consensus)}')

    return consensus_sequences[:count]


def main():

    if len(sys.argv) < 6:
        print("PARAMETRES: python main.py <vcf_file> <fasta_file> <LENGTH> <COUNT> <TRESHOLD> <output_file_name>")
        sys.exit(1)

    vcf_file = sys.argv[1]
    fasta_file = sys.argv[2]
    length = int(sys.argv[3])
    count = int(sys.argv[4])
    frequency_threshold = float(sys.argv[5])
    output_file = sys.argv[6]

    variants = parse_vcf(vcf_file)
    reference = parse_fasta(fasta_file)

    consensus_sequences = generate_consensus(
        variants, reference, length, count, frequency_threshold
    )

    with open(output_file, 'w') as out:
        out.write('\n'.join(consensus_sequences))

if __name__ == "__main__":
    main()