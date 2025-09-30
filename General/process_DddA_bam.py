import pysam
import argparse


"""
1) Correct aligned DddA BAM and replace likely deamination events with ambiguity codes (C|T: Y, G|A: R)
2) Create DA tag listing the moleular coordinates of deaminations
3) Add FD & LD tags for first and last deamination events in molecular coordinates
"""

# parse command line arguments
parser = argparse.ArgumentParser(description = "DddA BAM preprocessing",
    epilog = "")
parser.add_argument("-b", "--bam", required = True, metavar = '', help = "DddA aligned BAM to correct")
parser.add_argument("-c", "--sd_cutoff", required = False, default=3., type=float,  help = "Strand mut sigma probability cutoff [default:%(default)g]")
parser.add_argument(
        "-V",
        "--verbose",
        default=False,
        action="store_true",
        help="Be more verbose with output",
    )
args = parser.parse_args()


import logging

if args.verbose:
    logging.basicConfig(
            level=logging.INFO,
            format="%(message)s",
    )
    logging.info(str(args))


# identify fastq files in dir
bam_name = args.bam
sd_cutoff = args.sd_cutoff

assert sd_cutoff>0

def determine_da_strand_MD(read_obj, cutoff,do_logging=False):
    # based on the proportion of C->T & G->A determine the strand acted upon by DddA
    # only counting single base substitutions
    import math 
    seq = read_obj.query_sequence
    pair = read_obj.get_aligned_pairs(matches_only=False, with_seq=True)
    c = 0
    g = 0
    total = 0
    #for pos in pair:
    for qi,ri,ref_base in pair:
        if qi == None or ri == None: # indel, ignore
            pass
        else:
            ref_base = ref_base.upper()
            q_base = seq[qi]
            if q_base != ref_base:
                total += 1
                match (ref_base,q_base):
                    case ("C","T"):            
                        c += 1
                    case ("G","A"):
                        g += 1
                        
    
    
    cutoff_N = (c+g)/2 - cutoff*0.5*math.sqrt(c+g)  # 3 sigma below mean.
    #if c+g == 0:
    #    return('none')
    r="undetermined"
    if g<cutoff_N: #c/(c+g) >= cutoff:
        r=('CT')
    elif c<cutoff_N: # g/(c+g) >= cutoff:
        r=('GA')
    if do_logging:
        logging.info("%s:%d %s C:%d G:%d Cp:%g  %s cut:%g",read_obj.reference_name, read_obj.reference_start,read_obj.query_name,c,g,(0. if c+g==0 else c/(c+g)),r,cutoff_N)
    return r 
def check_num_assigned(sam_obj, cutoff):
    none = 0
    und = 0
    ct = 0
    ga = 0
    for read in sam_obj:
        if read.is_secondary == False and read.is_supplementary == False:
            change = determine_da_strand_MD(read, cutoff)
            if change == 'none':
                none += 1
            elif change == 'undetermined':
                und += 1
            elif change == 'CT':
                ct += 1
            else:
                ga += 1
    return({'CT':ct, 'GA':ga, 'Undetermined':und, 'None':none})

def correct_read_MD(read_obj, strand):
    """ Identify single-base changes from the reference that are likely induced by DddA
     and correct the original sequence using ambiguity codes (C|T: Y, G|A: R) and output new DA-tag positions.
     Limit detection to the previously identified DddA strand info (either ct or ga).
     Output everything in FIBER coordinates, not reference!
    """
    seq = read_obj.query_sequence
    pair = read_obj.get_aligned_pairs(matches_only=False, with_seq=True)
    new_seq = ''
    amb_codes = {'CT':'Y', 'GA':'R'}
    deam_pos = [] # mol coordinates of likely base changes
    for pos in pair:
        if pos[0] == None: # deletion, ignore
            pass
        elif pos[1] == None: # insertion, use seq base
            qi = pos[0]
            new_seq += seq[qi]
        else:
            qi = pos[0]
            ref_pos = pos[2].upper()
            if seq[qi] != ref_pos:
                change = ref_pos + seq[qi]
                if change == strand:
                    new_seq += amb_codes[change] # Update seq with ambiguity codes
                    deam_pos.append(qi+1) # track DA positions in 1-indexed
                else:
                    new_seq +=seq[qi]
            else:
                new_seq += seq[qi]
    return(new_seq, deam_pos)


# write corrected reads to new BAM
bam = pysam.AlignmentFile(bam_name, "rb",threads=4)

from pathlib import Path
new_bam = Path(bam_name).with_suffix(".corrected.bam").name
assert new_bam!=bam_name
corrected_bam = pysam.AlignmentFile(new_bam, "wb", template=bam,threads=4)
N_written=0
N_skipped=0
for read in bam:
    if read.is_secondary == False and read.is_supplementary == False and read.has_tag("MD"):
        strand = determine_da_strand_MD(read, sd_cutoff,do_logging=(N_written%10000)==0)
        MD = read.get_tag('MD')
        if strand in ['CT','GA']:
            # WRITE NEW seq with added tags
            new_seq, deam_pos = correct_read_MD(read, strand)
            quals = read.query_qualities
            read.query_sequence = new_seq
            read.seq = new_seq
            read.query_qualities = quals
            read.set_tags([('DA', deam_pos), ('FD', deam_pos[0], "i"), ('LD', deam_pos[-1], "i"), ('ST', strand), ('MD', MD)])
        else:
            read.set_tags([('DA', [0]), ('FD', 0, "i"), ('LD', 0, "i"), ('ST', strand),('MD', MD)])
        corrected_bam.write(read)
        N_written+=1
    else:
        N_skipped+=1
bam.close()
corrected_bam.close()
logging.info("Wrote %d reads, skipped %d (%g%%).",N_written,N_skipped,N_skipped/(N_skipped+N_written))
