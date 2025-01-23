import sys
import random
import time
import logging
from multiprocessing import Pool, cpu_count

def error_logger(func):
    """Декоратор для обработки ошибок."""
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception as e:
            logging.error(f"[ERROR]: {str(e)}")
            sys.exit(1)
    return wrapper

class VCFParser:
    """Класс для парсинга VCF-файлов."""
    def __init__(self, vcf_file):
        self.vcf_file = vcf_file
        self.variants = {}

    def parse(self):
        logging.info(f"Parsing VCF file: {self.vcf_file}")
        with open(self.vcf_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                cols = line.strip().split('\t')
                chrom, pos, ref, alt = cols[0], int(cols[1]), cols[3], cols[4]
                samples = cols[9:]
                if chrom not in self.variants:
                    self.variants[chrom] = []
                self.variants[chrom].append((pos, ref, alt, samples))
        for chrom in self.variants:
            self.variants[chrom].sort(key=lambda x: x[0])
        logging.info(f"Finished parsing and sorting VCF file: {self.vcf_file}")
        return self.variants


class FASTAParser:
    """Класс для парсинга FASTA-файлов."""
    def __init__(self, fasta_file):
        self.fasta_file = fasta_file
        self.sequences = {}

    def parse(self):
        logging.info(f"Parsing FASTA file: {self.fasta_file}")
        with open(self.fasta_file, 'r') as f:
            current_chrom = None
            for line in f:
                if line.startswith('>'):
                    current_chrom = line[1:].strip()
                    self.sequences[current_chrom] = []
                else:
                    self.sequences[current_chrom].append(line.strip())
        logging.info(f"Finished parsing FASTA file: {self.fasta_file}")
        return {chrom: ''.join(seq) for chrom, seq in self.sequences.items()}


class ConsensusGenerator:
    """Класс для генерации консенсусных последовательностей."""
    def __init__(self, variants, reference, length, count, frequency_threshold):
        self.variants = variants
        self.reference = reference
        self.length = length
        self.count = count
        self.frequency_threshold = frequency_threshold

    def generate_for_chromosome(self, chrom_data):
        chrom, sequence = chrom_data
        logging.info(f"Generating consensus for chromosome: {chrom}")
        if chrom not in self.variants:
            return []

        chrom_variants = self.variants[chrom]
        chrom_length = len(sequence)
        consensus_sequences = []

        current_count = self.count // len(self.reference)
        chrom_positions = set()

        while len(chrom_positions) < current_count:
            start = random.randint(0, chrom_length - self.length)
            if start in chrom_positions:
                continue
            chrom_positions.add(start)

            end = start + self.length
            consensus = list(sequence[start:end])
            for pos, ref, alt, samples in chrom_variants:
                if start <= pos < end:
                    sample_alleles = [s.split(':')[0] for s in samples if s]
                    total_samples = 2 * len(sample_alleles)
                    if sample_alleles:
                        allele_freq = sample_alleles.count('1') / total_samples
                        if allele_freq >= self.frequency_threshold:
                            consensus[pos - start] = alt

            consensus_sequences.append(f'>{chrom}_{start}_{end}\n{"".join(consensus)}')

        logging.info(f"Finished generating consensus for chromosome: {chrom}")
        return consensus_sequences

    def generate(self):
        logging.info("Generating consensus sequences using multiprocessing")
        with Pool(cpu_count()) as pool:
            results = pool.map(self.generate_for_chromosome, self.reference.items())
        consensus_sequences = [seq for sublist in results for seq in sublist]
        logging.info("Finished generating all consensus sequences")
        return consensus_sequences[:self.count]


@error_logger
def main():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    if len(sys.argv) < 7:
        logging.error("PARAMETERS: python main.py <vcf_file> <fasta_file> <LENGTH> <COUNT> <THRESHOLD> <output_file_name>")
        sys.exit(1)

    start_time = time.time()

    vcf_file, fasta_file, length, count, frequency_threshold, output_file = (
        sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5]), sys.argv[6]
    )

    vcf_parser = VCFParser(vcf_file)
    variants = vcf_parser.parse()
    fasta_parser = FASTAParser(fasta_file)
    reference = fasta_parser.parse()
    generator = ConsensusGenerator(variants, reference, length, count, frequency_threshold)
    consensus_sequences = generator.generate()

    with open(output_file, 'w') as out:
        out.write('\n'.join(consensus_sequences))

    end_time = time.time()
    logging.info(f"Execution Time: {end_time - start_time:.2f} seconds")

if __name__ == "__main__":
    main()