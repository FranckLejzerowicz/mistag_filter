import os
import subprocess
from multiprocessing import Process, Manager, current_process, cpu_count

from mistag_writer import write_fastas, write_merging_cmd
from mistag_utils import increment_nested

def get_merging_config(soft, folder):
    merging_config_fp = '%s/merging_%s.conf' % (folder, soft)
    merging_cmd_fp = '%s.cmd' % os.path.splitext(merging_config_fp)[0]
    if os.path.isfile(merging_cmd_fp):
        os.remove(merging_cmd_fp)
    merging = {}
    if os.path.isfile(merging_config_fp) == False:
        o=open(merging_config_fp, 'w')
        if soft == 'vsearch':
            o.write('min_qual,--fastq_qmin,0\n')
            o.write('max_qual,--fastq_qmax,41\n')
            o.write('min_overlap,--fastq_minovlen,15\n')
            o.write('max_diffs,--fastq_maxdiffs,5\n')
            o.write('min_size,--fastq_minmergelen,1\n')
            o.write('max_size,--fastq_maxmergelen,1000\n')
            o.write('threads,--threads,%s\n' % cpu_count())
        elif soft == 'pandaseq':
            o.write('thresh,-t,0.6\n')
            o.write('algo,-A,pear\n')
            o.write('threads,-T,%s\n' % cpu_count())
            o.write('min_len,-l,\n')
            o.write('max_len,-L,\n')
            o.write('min_overap,-o,30\n')
            o.write('comment,-d,bfsrk\n')
        o.close()
    if soft == 'vsearch':
        merging['for'] = ['--fastq_mergepairs', None]
        merging['rev'] = ['--reverse', None]
        merging['out'] = ['--fastaout', None]
    elif soft == 'pandaseq':
        merging['for'] = ['-f', None]
        merging['rev'] = ['-r', None]
        merging['out'] = ['-w', None]
        merging['log'] = ['-G', None]
    with open(merging_config_fp) as f:
        for line in f:
            line_split = line.strip().split(',')
            merging[line_split[0]] = line_split[1:]
    return merging, merging_cmd_fp


def run_merging(folder, fastqs, merging_soft, merging_config, return_dict,
                merging_cmd_fp):
    curProc_name = current_process().name
    cmd = [merging_soft]
    for config_opt, config_val in merging_config.items():
        if config_val[-1]:
            cmd.extend(config_val)
        else:
            if config_opt == 'for':
                cmd.extend([config_val[0], fastqs[0]])
            if config_opt == 'rev':
                cmd.extend([config_val[0], fastqs[1]])
            if config_opt == 'log':
                cmd.extend([config_val[0], '%s_log.bz2' %
                            os.path.splitext(fastqs[0])[0]])
            if config_opt == 'out':
                out = '%s.fasta' % '_'.join(fastqs[0].split('_')[:-1])
                return_dict[out] = [curProc_name]
                cmd.extend([config_val[0], out])
    process = subprocess.Popen(cmd)
    process.wait()
    write_merging_cmd(merging_cmd_fp, cmd)


def get_merging_stats(return_dict, folder, stats_merging):
    """Collect the number of reads in each the fastqs
    sequences submitted to the merging process
    """
    grep_key = ''
    # return_dict: {'<sample>.fasta': 'merging_<sample>', ...}
    # get grep_key as common sequencing identifier (i.e. flowcell)
    for filin in return_dict.keys():
        with open(filin) as f:
            for line in f:
                grep_key = '\:'.join(line[1:].split(':')[:3])
                break
        if grep_key:
            break
    # grep this each-sequence-entry identifier in each demultiplexed fastq
    cmd = 'grep -cR "%s" %s/*_fwd.fastq' % (grep_key, folder)
#    try:
#        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, encoding='utf8')
#    except TypeError:
#        p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    p.wait()
    # get for each sample the number of sequences submitted to merging
    for fastq_grep in p.stdout:
        fastq_grep_decode = fastq_grep.decode('utf8').strip()
        f = str(fastq_grep_decode).strip().split(':')[0]
        n = str(fastq_grep_decode).strip().split(':')[1]
        file_key = os.path.abspath('%s.fasta' % f.split('_fwd.fastq')[0])
        stats_merging[file_key] = [n]


def get_derep(return_dict, samples_fastqs, mistag_unexpected, stats_merging):
    derep = {}
    # for each fasta file outputed by the merging algorithm
    for fasta in return_dict.keys():
        # get the relevant deign combination (except if mistag)
        if fasta.endswith('mistag.fasta'):
            combi = 0
        else:
            combi = samples_fastqs[fasta.replace('.fasta', '')][-1]
        # parse the fasta file
        c = 0
        add = 0
        seq = ''
        with open(fasta) as f:
            for line in f:
                # if header
                if line[0] == '>':
                    # once previous sequence consumed
                    if len(seq):
                        # add previous sequence combi info if told to be added
                        if add:
                            c += 1
                            increment_nested(derep, seq, combi)
                        seq = ''
                    # every sequence is to be added
                    add = 1
                    # but not necessarily for mistags
                    if fasta.endswith('_mistag.fasta'):
                        # add only if sequence is in a plausible sample
                        ID = ':'.join(line[1:].split(':')[:7])
                        if ID in mistag_unexpected:
                            combi = mistag_unexpected[ID]
                        else:
                            add = 0
                # collect sequence
                else:
                    seq += line.strip()
            # process the last sequence entry of the fasta file
            if add:
                c += 1
                increment_nested(derep, seq, combi)
            # update stats_count if the fasta file contain the merged mistags
            if fasta.endswith('_mistag.fasta') == False:
                file_key = os.path.abspath(fasta)
                stats_merging[file_key].append(str(c))
    return derep


def perform_merging(args, fastin, samples_fastqs, mistag_fastqs,
                    stats_merging, mistag_unexpected, multiproc):
    merging_soft = args['m']
    merging_config, merging_cmd_fp = get_merging_config(merging_soft, fastin)
    if multiproc:
        jobs_merging = []
        manager = Manager()
        return_dict = manager.dict()
        # perform the merging for each sample fastqs
        for sample_fastq, sample_fastqs in samples_fastqs.items():
            pname = 'merging_%s' % sample_fastq.split('/')[-1]
            print('Process:', pname)
            p = Process(target = run_merging, args = (fastin,
                                                    sample_fastqs,
                                                    merging_soft,
                                                    merging_config,
                                                    return_dict,
                                                    merging_cmd_fp,),
                                                    name=pname)
            jobs_merging.append(p)
            p.start()
        for job_merging in jobs_merging:
            job_merging.join()
    else:
        return_dict = {}
        for sample_fastq, sample_fastqs in samples_fastqs.items():
            pname = 'merging_%s' % sample_fastq.split('/')[-1]
            print('Process:', pname)
            run_merging(fastin, sample_fastqs, merging_soft, merging_config,
                        return_dict, merging_cmd_fp)
    # also perform the merging for the mistag fastqs
    pname = 'merging_mistags'
    print('Process:', pname)
    run_merging(fastin, mistag_fastqs, merging_soft, merging_config,
                return_dict, merging_cmd_fp)
    print('All mergings complete!')
    get_merging_stats(return_dict, fastin, stats_merging)
    derep = get_derep(return_dict, samples_fastqs, mistag_unexpected,
                      stats_merging)
    fastout = write_fastas(fastin, derep, stats_merging)
    return fastout

