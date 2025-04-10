#!/usr/bin/env python
# conda create -c bioconda -n pysam_env pysam
# For BAM/SAM files: python script.py input.bam -o junctions.tsv
# For CRAM files (with reference): python script.py input.cram -r reference.fa -o junctions.tsv
# With custom mapping quality threshold: python script.py input.bam -q 30 -o junctions.tsv
# Output to stdout: python script.py input.bam


import pysam
import sys
import re   
import argparse

def main():
    parser = argparse.ArgumentParser(description='Process SAM/BAM/CRAM to detect split reads and count junctions')
    parser.add_argument('input_file', help='Input SAM/BAM/CRAM file')
    parser.add_argument('-q', '--min-mapq', type=int, default=255,
                       help='Minimum mapping quality (default: 255)')
    parser.add_argument('-r', '--reference', help='Reference genome FASTA (required for CRAM)')
    parser.add_argument('-o', '--output', help='Output file (default: stdout)')
    args = parser.parse_args()

    HoS = {}
    minmapq = args.min_mapq

    try:
        # Determine file mode and open appropriately
        file_mode = "rb"  # Binary mode for BAM/CRAM
        open_kwargs = {}
        
        if args.input_file.endswith('.cram'):
            if not args.reference:
                print("Error: Reference genome required for CRAM files", file=sys.stderr)
                sys.exit(1)
            open_kwargs['reference_filename'] = args.reference
            file_mode = "rc"  # CRAM specific mode
        
        with pysam.AlignmentFile(args.input_file, file_mode, **open_kwargs) as fh:
            for read in fh:
                if read.mapping_quality >= minmapq and read.cigarstring and 'N' in read.cigarstring:
                    chr = read.reference_name
                    start = read.reference_start + 1  # Convert to 1-based
                    cigar = read.cigarstring
                    
                    # Parse CIGAR string
                    inter = [int(x) for x in re.split(r'\D', cigar) if x]
                    cops = [x for x in re.split(r'\d', cigar) if x]
                    
                    for i in range(len(cops)):
                        ofs = 0
                        if cops[i] == 'N':
                            for k in range(i):
                                if cops[k] in ['M', 'N', 'D']:
                                    ofs += inter[k]
                                elif cops[k] == 'I':
                                    ofs -= inter[k]
                                # Ignore S and s
                            
                            junc_start = start + ofs
                            junc_end = start + ofs + inter[i] - 1
                            
                            hk = f"{chr}_{junc_start}_{junc_end}"
                            HoS[hk] = HoS.get(hk, 0) + 1

    except Exception as e:
        print(f"Error processing file: {e}", file=sys.stderr)
        sys.exit(1)

    # Output results
    output_handle = open(args.output, 'w') if args.output else sys.stdout
    try:
        for key, counts in HoS.items():
            output_handle.write("\t".join(key.split('_') + [str(counts)]) + "\n")
    finally:
        if args.output:
            output_handle.close()

if __name__ == "__main__":
    main()