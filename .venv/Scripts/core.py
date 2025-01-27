import random
import logging
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
        self.chrom_lengths = {}
        self.variants = {}

    @log_exceptions
    def parse_reference_lengths(self):
        with open(self.reference_path, 'r') as file:
            current_chrom, current_length = None, 0
            for line in file:
                if line.startswith('>'):
                    if current_chrom:
                        self.chrom_lengths[current_chrom] = current_length
                    current_chrom = line.strip()[1:]
                    current_length = 0
                else:
                    current_length += len(line.strip())
            if current_chrom:
                self.chrom_lengths[current_chrom] = current_length
        logger.info(f"Parsed chromosome lengths: {self.chrom_lengths}")

    @log_exceptions
    def parse_vcf(self):
        with open(self.vcf_path, 'r') as file:
            sample_index = None
            for line in file:
                if line.startswith("#CHROM"):
                    headers = line.strip().split('\t')
                    if self.sample_name in headers:
                        sample_index = headers.index(self.sample_name)
                        logger.info(f"Found sample '{self.sample_name}' at index {sample_index}")
                    else:
                        raise ValueError(f"Sample '{self.sample_name}' not found in VCF headers.")
                elif not line.startswith('#') and sample_index is not None:
                    fields = line.strip().split('\t')
                    chrom, pos, ref, alt = fields[0], int(fields[1]) - 1, fields[3], fields[4]
                    genotype = fields[sample_index].split(':')[0]
                    if "1" in genotype:
                        self.variants.setdefault(chrom, []).append((pos, ref, alt))
        if sample_index is None:
            raise ValueError(f"Sample '{self.sample_name}' not found in the VCF file.")
        logger.info(f"Parsed variants for sample {self.sample_name}")

    @log_exceptions
    def generate_start_positions(self):
        all_positions = []
        for chrom, length in self.chrom_lengths.items():
            max_start = max(0, length - self.length_consensus + 1)
            if max_start <= 0:
                logger.warning(f"Consensus length {self.length_consensus} exceeds chromosome {chrom} length {length}. Skipping...")
                continue

            if self.random_start:
                positions = random.sample(range(max_start), min(self.total_consensus, max_start))
            else:
                positions = list(range(0, max_start, self.length_consensus))[:self.total_consensus]
            all_positions.extend([(chrom, pos) for pos in positions])

        if len(all_positions) > self.total_consensus:
            all_positions = random.sample(all_positions, self.total_consensus)
        logger.info(f"Generated start positions: {all_positions}")
        return all_positions

    @log_exceptions
    def create_consensus_sequences(self, start_positions):
        consensus_sequences = {}
        chrom_sequences = self._load_reference_sequences()

        for chrom, start in start_positions:
            original_seq = chrom_sequences[chrom][start:start + self.length_consensus]
            modified_seq = self.apply_variants_to_sequence(chrom, start, original_seq)
            change_label = "MODIFIED" if original_seq != modified_seq else "UNMODIFIED"
            header = f"{self.sample_name}_{chrom}_Chunk_{start + 1}-{start + self.length_consensus}_{change_label}"
            consensus_sequences[header] = modified_seq

        logger.info(f"Created {len(consensus_sequences)} consensus sequences")
        return consensus_sequences

    def _load_reference_sequences(self):
        sequences = {}
        with open(self.reference_path, 'r') as file:
            current_chrom, current_seq = None, []
            for line in file:
                if line.startswith('>'):
                    if current_chrom:
                        sequences[current_chrom] = ''.join(current_seq)
                    current_chrom = line.strip()[1:]
                    current_seq = []
                else:
                    current_seq.append(line.strip())
            if current_chrom:
                sequences[current_chrom] = ''.join(current_seq)
        return sequences

    def apply_variants_to_sequence(self, chrom, start, sequence):
        if chrom not in self.variants:
            return sequence

        seq = list(sequence)
        for pos, ref, alt in self.variants[chrom]:
            if start <= pos < start + len(sequence):
                local_pos = pos - start
                if seq[local_pos:local_pos + len(ref)] == list(ref):
                    seq[local_pos:local_pos + len(ref)] = list(alt)
                else:
                    logger.warning(f"Variant mismatch at {chrom}:{pos + 1}. Skipping variant.")
        return ''.join(seq)

    @log_exceptions
    def write_fasta(self, sequences):
        with open(self.output_path, 'w') as fasta:
            for header, seq in sequences.items():
                fasta.write(f">{header}\n{seq}\n")
        logger.info(f"FASTA written to {self.output_path}")

    @log_exceptions
    def run(self):
        start_time = time.time()
        self.parse_reference_lengths()
        self.parse_vcf()
        start_positions = self.generate_start_positions()
        consensus_sequences = self.create_consensus_sequences(start_positions)
        self.write_fasta(consensus_sequences)
        logger.info(f"Total execution time: {time.time() - start_time:.2f} seconds")

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

    generator = ConsensusGenerator(
        args.reference, args.vcf, args.output, args.length, args.total, args.random, args.sample
    )
    generator.run()

