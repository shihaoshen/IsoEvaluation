import sys, os, re
import gffutils as gu
from collections import Counter
from intervaltree import Interval, IntervalTree

def is_cons(tid, anno_trans, iso_trans, head_constr, tail_constr):
    """
    :param trans1: [(exon1), (exon2), (exon3)]
    :param trans2: [(exon1), (exon2), (exon3)]
    :return: YES(consistent), NO(in-consistent)
    """

    anno_start = anno_trans[0][0]
    anno_end = anno_trans[-1][1]
    iso_start = iso_trans[0][0]
    iso_end = iso_trans[-1][1]

    if __debug__:
        print 'before:\n', tid, anno_trans, iso_trans, iso_start, iso_end, head_constr, tail_constr

    if anno_start > iso_start or anno_end < iso_end:
        return False

    anno = IntervalTree(Interval(*iv) for iv in anno_trans)
    iso = IntervalTree(Interval(*iv) for iv in iso_trans)

    anno.merge_overlaps()
    iso.merge_overlaps()

    if not head_constr: anno.chop(0, iso_start)
    if not tail_constr: anno.chop(iso_end, anno_end)

    if __debug__:
        print 'after:\n', tid, anno, iso
        if anno == iso:
            print 'match:\n', tid, iso_start, iso_end, anno, iso

    return anno == iso


def cons_trans_from_gtf(gtf_db, trans_abund, chr, exon_coor, isoform, iso_i):
    abund = 0

    #gene = gtf_db[gene_name]
    iso_coor = []

    if __debug__:
        print 'iso_i:', iso_i

    head_constr = (isoform[1] != 0)
    tail_constr = (isoform[-1] != len(exon_coor)-1)
    for i in range(isoform[0]):
        iso_coor.append((exon_coor[isoform[i + 1]][0] - 1, exon_coor[isoform[i + 1]][1]))

    gname = []
    if __debug__:
        print 'region:', chr, exon_coor[0][0], exon_coor[-1][1]
    for g in gtf_db.region((chr, exon_coor[0][0], exon_coor[-1][1]), completely_within=False, featuretype=['gene']):
        gname.append(g.attributes['gene_id'][0])

    for g in gname:
        gene = gtf_db[g]
        for trans in gtf_db.children(gene, featuretype='transcript', order_by='start'):
            trans_coor = []
            # print trans.attributes['transcript_id'][0]
            for exon in gtf_db.children(trans, featuretype='exon', order_by='start'):
                trans_coor.append((exon.start - 1, exon.end))
            if is_cons(trans.attributes['transcript_id'][0], trans_coor, iso_coor, head_constr, tail_constr):
                # print 'abund', trans.attributes['transcript_id'][0], trans_abund[trans.attributes['transcript_id'][0]]
                abund += trans_abund[trans.attributes['transcript_id'][0]]

    return abund


def print_abund(asm_name, iso_abund):
    if iso_abund:
        print asm_name, len(iso_abund),

        for a in iso_abund:
            if sum(iso_abund) != 0:
                print float(a) / sum(iso_abund),
            else:
                print 0.0,

        for a in iso_abund:
            print a,

        print

def gen_psi_from_iso_fp(iso_fp, trans_abund, gtf_db):
    line = iso_fp.readline()
    #gene_name = ''
    exon_coor = []
    iso_abund = []
    isoform = []
    asm_name = ''
    chr = ''
    while line:
        if line[0] == 'A':
            # 5. calculate accumulative abundance for each isoform
            if isoform:
                iso_i = 0
                for iso in isoform:
                    iso_i += 1
                    # 4. for each isoform, sum up all the consistent trans
                    iso_abund.append(cons_trans_from_gtf(gtf_db, trans_abund, chr, exon_coor, iso, iso_i))

            # 5. output relative ratio for each isoform
            print_abund(asm_name, iso_abund)
            iso_abund = []
            isoform = []

            asm_name = line.split()[0]
            chr = line.split()[4]
            #gene_name = line.split()[6]
            exon_line = iso_fp.readline().split()
            exon_coor = []
            for e in exon_line:
                exon_coor.append((int(e.split(',')[0]), int(e.split(',')[1])))

        else:  # isoform structure
            isoform.append(map(int, line.split()))

        line = iso_fp.readline()

    if isoform:
        iso_i = 0
        for iso in isoform:
            iso_i += 1
            # 4. for each isoform, sum up all the consistent trans
            iso_abund.append(cons_trans_from_gtf(gtf_db, trans_abund, chr, exon_coor, iso, iso_i))

    print_abund(asm_name, iso_abund)


def gen_psi_from_rmats_fp(rmats_fp, trans_abund, gtf_db):
    firstline = rmats_fp.readline().split()
    HEADER = {firstline[i]: i for i in range(len(firstline))}

    # print HEADER
    for line in rmats_fp:
        exon_coor = []
        iso_abund = []
        ele = line.split()
        asm_name = ele[HEADER['ID']]
        gene_name = ele[HEADER['Gene_ID']]
        exon_coor = [
            (int(ele[HEADER['Up_start']]), int(ele[HEADER['Up_end']])),
            (int(ele[HEADER['SE_start']]), int(ele[HEADER['SE_end']])),
            (int(ele[HEADER['Down_start']]), int(ele[HEADER['Down_end']])),
        ]
        # 5. calculate accumulative abundance for each isoform
        isoform = [3, 0, 1, 2]
        iso_i = 0
        iso_abund.append(cons_trans_from_gtf(gtf_db, trans_abund, gene_name, exon_coor, isoform, iso_i))
        isoform = [2, 0, 2]
        iso_abund.append(cons_trans_from_gtf(gtf_db, trans_abund, gene_name, exon_coor, isoform, iso_i))

        # 6. output relative ratio for each isoform
        if iso_abund:
            print asm_name, len(iso_abund),

            for a in iso_abund:
                if sum(iso_abund) != 0:
                    print float(a) / sum(iso_abund),
                else:
                    print 0.0,

            for a in iso_abund:
                print a,
            print
            # print ele[HEADER['PS1']], 1- ele[HEADER['PS1']], ele[HEADER['PS2']], 1- ele[HEADER['PS2']]


if __name__ == '__main__':
    if len(sys.argv) != 4:
        print "Usage:"
        print "%s in.pro in.GTF in.IsoExon > psi.out" % sys.argv[0]
        print "%s in.pro in.GTF in.JC.psi > psi.out" % sys.argv[0]
        sys.exit(1)

    pro_fp = open(sys.argv[1])
    gtf_fn = sys.argv[2]
    iso_fp = None
    rmats_fp = None
    if re.search(r'Exon$', sys.argv[3]):
        iso_fp = open(sys.argv[3])
    elif re.search(r'psi$', sys.argv[3]):
        rmats_fp = open(sys.argv[3])

    # 1. read Flux profile file
    #print('read Flux profile file');
    trans_abund = Counter()
    line = pro_fp.readline()
    while line:
        trans_abund[line.split()[1]] = int(line.split()[5])
        line = pro_fp.readline()

    # 2. read gene_name and exon_coor for each isoform
    #print('rread gene_name and exon_coor for each isoform');
    if os.path.isfile(gtf_fn + '.gffdb'):
        gtf_db = gu.FeatureDB(gtf_fn + '.gffdb')
    else:
        gtf_db = gu.create_db(gtf_fn, gtf_fn + '.gffdb', disable_infer_genes = True, disable_infer_transcripts = True)

    # 3. obtain consistent transcripts
    #print('obtain consistent transcripts')
    if iso_fp:
        gen_psi_from_iso_fp(iso_fp, trans_abund, gtf_db)
    elif rmats_fp:
        gen_psi_from_rmats_fp(rmats_fp, trans_abund, gtf_db)

    pro_fp.close()
    if iso_fp: iso_fp.close()
    if rmats_fp: rmats_fp.close()
