import sys
import random

def error_logger(func):
    """Декоратор для вывода ошибок"""
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception as e:
            print("\n[ERROR]:", str(e))
            sys.exit(1)
    return wrapper

def parse_vcf(vcf_file):
    """Парсит VCF-файл и возвращает словарь с вариантами"""
    variants = {}
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            cols = line.strip().split('\t')
            chrom, pos, ref, alt = cols[0], int(cols[1]), cols[3], cols[4]
            samples = cols[9:]
            variants.setdefault(chrom, []).append((pos, ref, alt, samples))
    return variants

def parse_fasta(fasta_file):
    """Парсит FASTA-файл и возвращает словарь с референсным геномом"""
    sequences = {}
    with open(fasta_file, 'r') as f:
        current_chrom = None
        for line in f:
            if line.startswith('>'):
                current_chrom = line[1:].strip()
                sequences[current_chrom] = []
            else:
                sequences[current_chrom].append(line.strip())
    return {chrom: ''.join(seq) for chrom, seq in sequences.items()}

def generate_consensus(variants, reference, length, count, frequency_threshold):
    """Генерирует консенсусные последовательности."""
    consensus_sequences = []
    chrom_count = len(reference)
    per_chrom_count, extra = divmod(count, chrom_count)

    for i, (chrom, sequence) in enumerate(reference.items()):
        if chrom not in variants:
            continue

        chrom_variants = variants[chrom]
        chrom_length = len(sequence)

        if length > chrom_length:
            raise ValueError(f"Length of consensus ({length}) is greater than the chromosome length ({chrom_length})")

        current_count = per_chrom_count + (1 if i < extra else 0)
        chrom_positions = set()

        while len(chrom_positions) < current_count:
            start = random.randint(0, chrom_length - length)
            if start in chrom_positions:
                continue
            chrom_positions.add(start)

            end = start + length
            consensus = list(sequence[start:end])
            for pos, ref, alt, samples in chrom_variants:
                if start <= pos < end:
                    sample_alleles = [s.split(':')[0] for s in samples if s]
                    if sample_alleles:
                        allele_freq = sample_alleles.count('1') / len(sample_alleles)
                        if allele_freq >= frequency_threshold:
                            consensus[pos - start] = alt

            consensus_sequences.append(f'>{chrom}_{start}_{end}\n{"".join(consensus)}')

    return consensus_sequences[:count]

@error_logger
def main():
    if len(sys.argv) < 7:
        print("PARAMETERS: python main.py <vcf_file> <fasta_file> <LENGTH> <COUNT> <THRESHOLD> <output_file_name>")
        sys.exit(1)

    vcf_file, fasta_file, length, count, frequency_threshold, output_file = (
        sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5]), sys.argv[6]
    )

    variants = parse_vcf(vcf_file)
    reference = parse_fasta(fasta_file)
    consensus_sequences = generate_consensus(variants, reference, length, count, frequency_threshold)

    with open(output_file, 'w') as out:
        out.write('\n'.join(consensus_sequences))

main()