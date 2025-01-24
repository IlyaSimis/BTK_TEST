import random
import logging
import multiprocessing
import argparse
import time
from functools import wraps

logging.basicConfig(level=logging.INFO, format='%(asctime)s [%(levelname)s] %(message)s')
logger = logging.getLogger(__name__)

def log_exceptions(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except Exception as e:
            logger.exception(f"Error in function {func.__name__}")
            raise e
    return wrapper

class ConsensusGenerator:
    def __init__(self, reference_path, vcf_path, output_path, length_consensus, total_consensus, random_start, sample_name):
        self.reference_path = reference_path
        self.vcf_path = vcf_path
        self.output_path = output_path
        self.length_consensus = length_consensus
        self.total_consensus = total_consensus
        self.random_start = random_start
        self.sample_name = sample_name
        self.reference = None
        self.modified_reference = None
        self.variants = None

    @log_exceptions
    def load_fasta(self):
        reference = {}
        current_chrom = None
        with open(self.reference_path, 'r') as file:
            for line in file:
                if line.startswith('>'):
                    current_chrom = line.strip()[1:]
                    reference[current_chrom] = []
                elif current_chrom:
                    reference[current_chrom].append(line.strip())
        self.reference = {chrom: ''.join(seq) for chrom, seq in reference.items()}
        for chrom, seq in self.reference.items():
            logger.info(f"Loaded chromosome {chrom}: {len(seq)} bp")

    @log_exceptions
    def parse_vcf(self):
        variants = {}
        sample_index = None
        with open(self.vcf_path, 'r') as file:
            for line in file:
                if line.startswith("#CHROM"):
                    headers = line.strip().split('\t')
                    if self.sample_name in headers:
                        sample_index = headers.index(self.sample_name)
                    else:
                        raise ValueError(f"Sample '{self.sample_name}' not found in VCF headers.")
                elif not line.startswith('#'):
                    fields = line.strip().split('\t')
                    chrom = fields[0]
                    pos = int(fields[1]) - 1
                    ref = fields[3]
                    alt = fields[4]
                    if sample_index is not None:
                        genotype = fields[sample_index].split(':')[0]
                        if "1" in genotype:
                            if chrom not in variants:
                                variants[chrom] = []
                            variants[chrom].append((pos, ref, alt))
        self.variants = variants
        for chrom, vars in self.variants.items():
            logger.info(f"Parsed {len(vars)} variants for chromosome {chrom} for sample {self.sample_name}")

    @log_exceptions
    def apply_variants(self):
        self.modified_reference = self.reference.copy()
        for chrom, vars in self.variants.items():
            seq = list(self.modified_reference[chrom])
            for pos, ref, alt in vars:
                if seq[pos:pos + len(ref)] == list(ref):
                    seq[pos:pos + len(ref)] = list(alt)
            self.modified_reference[chrom] = ''.join(seq)
            logger.info(f"Variants applied to chromosome {chrom}")

    @log_exceptions
    def generate_start_positions(self):
        chrom_lengths = {chrom: len(seq) for chrom, seq in self.modified_reference.items()}
        all_positions = []
        for chrom, length in chrom_lengths.items():
            max_start = length - self.length_consensus + 1
            if max_start <= 0:
                logger.warning(f"Consensus length {self.length_consensus} exceeds chromosome {chrom} length {length}. Skipping...")
                continue
            if self.random_start:
                chrom_positions = random.sample(range(0, max_start), min(self.total_consensus, max_start))
            else:
                chrom_positions = list(range(0, max_start, self.length_consensus))[:self.total_consensus]
            all_positions.extend([(chrom, pos) for pos in chrom_positions])
        if len(all_positions) > self.total_consensus:
            all_positions = random.sample(all_positions, self.total_consensus)
        logger.info(f"Generated start positions: {all_positions}")
        return all_positions

    @log_exceptions
    def create_consensus_sequences(self, start_positions):
        consensus_sequences = {}
        for i, (chrom, start) in enumerate(start_positions):
            original_seq = self.reference[chrom][start:start + self.length_consensus]
            modified_seq = self.modified_reference[chrom][start:start + self.length_consensus]
            change_label = "MODIFIED" if original_seq != modified_seq else "UNMODIFIED"
            start_seq = start + 1
            end_seq = start + self.length_consensus
            header = f"{self.sample_name}_{chrom}_Chunk_{i+1}_({start_seq}-{end_seq})_{change_label}"
            consensus_sequences[header] = modified_seq
            logger.info(f"Created consensus sequence for {header}: {modified_seq}")
        return consensus_sequences

    @log_exceptions
    def write_fasta(self, sequences):
        with open(self.output_path, 'w') as fasta:
            for header, seq in sequences.items():
                fasta.write(f">{header}\n")
                fasta.write(f"{seq}\n")
        logger.info(f"FASTA written to {self.output_path}")

    @log_exceptions
    def run(self):
        start_time = time.time()
        self.load_fasta()
        self.parse_vcf()
        self.apply_variants()
        start_positions = self.generate_start_positions()
        consensus_sequences = self.create_consensus_sequences(start_positions)
        self.write_fasta(consensus_sequences)
        elapsed_time = time.time() - start_time
        logger.info(f"Total execution time: {elapsed_time:.2f} seconds")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate consensus FASTA from VCF and reference genome.")
    parser.add_argument("--reference", required=True, help="Path to the reference FASTA file.")
    parser.add_argument("--vcf", required=True, help="Path to the VCF file.")
    parser.add_argument("--output", required=True, help="Path to the output FASTA file.")
    parser.add_argument("--length", type=int, required=True, help="Length of consensus sequences.")
    parser.add_argument("--total", type=int, required=True, help="Total number of consensus sequences.")
    parser.add_argument("--random", action="store_true", help="Generate start positions randomly.")
    parser.add_argument("--sample", required=True, help="Sample name to include in FASTA headers.")

    args = parser.parse_args()

    process = multiprocessing.Process(target=ConsensusGenerator(
        args.reference, args.vcf, args.output, args.length, args.total, args.random, args.sample
    ).run)
    process.start()
    process.join()