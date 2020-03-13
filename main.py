import argparse
import os.path
import subprocess
from multiprocessing import *
import psutil

cpumemuse=0.3


def DNA_3GS(args):
    Ref=args.argR
    FastqFolder=args.argF
    SoftwareSelector=args.argS
    print("DNA_3GS",Ref,FastqFolder,SoftwareSelector)
    if SoftwareSelector == 32:
        print("Error, PBMM2 with MD tag can not recognized by sniffles")
        return
    #SoftwareSelector
    # tens --- map 1:minmap2    2:NGMLR     3:PBMM2
    # units --- variant 1:sniffles      2:PBSV

    Refpath, tempfilename = os.path.split(Ref);
    #shotname, extension = os.path.splitext(tempfilename);
    CPUcount = cpu_count()
    mem = psutil.virtual_memory()
    TotalMem = int(mem.total) / 1024 / 1024 / 1024
    print(Refpath,CPUcount,TotalMem)

    if SoftwareSelector // 10 == 1:  #minimap2
        #create ref.mmi
        # minimap2 -d ref.a.mmi human_g1k_v37.fast
        if not os.path.isfile(os.path.join(Refpath,'ref.mmi')):
            cmdlist = ['minimap2','-d',os.path.join(Refpath,'ref.mmi'),Ref]
            if subprocess.check_call(cmdlist)!=0:
                print("ref.mmi is fail")
                return
        print("ref.mmi is OK!")
        # fastq---->sam
        # minimap2 -MD -ax map-pb ref.fa read1.fa read2.fa read3.fa > merged.sam
        cmdlist = ['minimap2', '--MD','--cs', '-ax', 'map-pb', os.path.join(Refpath,'ref.mmi')]
        files = os.listdir(FastqFolder)
        for fqfile in files:
            if fqfile.endswith('fastq') or fqfile.endswith('fq'):
                cmdlist.extend([os.path.join(FastqFolder,fqfile)])
        #cmdlist.extend([os.path.join(FastqFolder, "ERRcp.fq")])
        cmdlist.extend(['>', 'merged.sam'])
        if subprocess.call(' '.join(cmdlist), shell=True)!= 0:
            print("minimap2: map is fail")
            return
        print("minimap2: map is OK!")
        # samtools sam---> sorted.bam
        cmdlist = ['samtools', 'view', '-bS', 'merged.sam', '>', 'merged.bam']
        if subprocess.call(' '.join(cmdlist), shell=True) != 0:
            print("samtools: convert is fail")
            return
        cmdlist = ['samtools', 'sort', '-@', str(max(int(CPUcount*cpumemuse),4)), '-m', str(max(int(TotalMem*cpumemuse),4))+'G', '-O', 'bam', '-o','merged.sorted.bam','merged.bam']
        if subprocess.call(' '.join(cmdlist), shell=True) != 0:
            print("samtools: sorted is fail")
            return
        print("samtools: sorted is OK!")
    elif SoftwareSelector // 10 == 2:    #NGMLR
        #ngmlr -t 10 -x pacbio -r /hadoop/quz/DNA_3GS/ref_data_cp/human_g1k_v37.fasta.gz -q /hadoop/quz/DNA_3GS/fastq_data_cp/ERRcp.fq
        files = os.listdir(FastqFolder)
        fqcount=1
        for fqfile in files:
            if fqfile.endswith('fastq') or fqfile.endswith('fq'):
                cmdlist = ['ngmlr', '-t',str(max(int(CPUcount*cpumemuse),4)), '-x', 'pacbio',
                           '-r', Ref, '-q', os.path.join(FastqFolder,fqfile), '-o',os.path.join(FastqFolder,'sub.{}.sam'.format(fqcount))]
                print("ngmlr --------- list: ",cmdlist)
                if subprocess.call(' '.join(cmdlist), shell=True) != 0:
                    print("NGMLR: map is fail")
                    return
            fqcount+=1
        print("NGMLR: map is OK!")
        # samtools sam---> sorted.bam
        files = os.listdir(FastqFolder)
        for fqfile in files:
            if fqfile.endswith('sam'):
                shotname, extension = os.path.splitext(fqfile);
                cmdlist = ['samtools', 'view', '-bS', os.path.join(FastqFolder,fqfile), '>', os.path.join(FastqFolder,shotname+'.bam')]
                if subprocess.call(' '.join(cmdlist), shell=True) != 0:
                    print("samtools: convert is fail")
                    return
        # merge all bam
        cmdlist=['samtools','merge','merged.sorted.bam']
        files = os.listdir(FastqFolder)
        for fqfile in files:
            if fqfile.endswith('bam'):
                cmdlist.extend([os.path.join(FastqFolder,fqfile)])
        if subprocess.call(' '.join(cmdlist), shell=True) != 0:
            print("samtools: merge is fail")
            return
        print("samtools: merge is OK!")

    elif SoftwareSelector // 10 == 3:    #PBMM2
        # create ref.mmi
        # pbmm2 index human_g1k_v37.fast ref.a.mmi
        if not os.path.isfile(os.path.join(Refpath, 'ref.mmi')):
            cmdlist = ['pbmm2', 'index', Ref, os.path.join(Refpath, 'ref.mmi')]
            if subprocess.check_call(cmdlist) != 0:
                print("ref.mmi is fail")
                return
        print("ref.mmi is OK!")

        # fastq ---> merge.sort.bam
        #create tmp file with multiple fq file path
        with open('tmplistfq.fofn','w') as fw:
            files = os.listdir(FastqFolder)
            for fqfile in files:
                if fqfile.endswith('fastq') or fqfile.endswith('fq'):
                    fw.write(os.path.join(FastqFolder,fqfile)+'\n')

        cmdlist = ['pbmm2','align',os.path.join(Refpath,'ref.mmi'),'tmplistfq.fofn','merged.sorted.bam','-j',str(max(int(CPUcount*cpumemuse),4)),
                   '-J',str(max(int(CPUcount*cpumemuse),4)),'--sort', '--preset', 'CCS', '--median-filter']
        if subprocess.call(' '.join(cmdlist), shell=True)!= 0:
            print("pbmm2 map is fail")
            return
        print("pbmm2 map is OK!")

    ###################################################################
    if SoftwareSelector % 10 == 1:  # sniffles
        # sniffles sort.bam--->vcf
        # sniffles -s 3 -z 2 -t 4 --report_seq --genotype -m aligned.sort.bam -v out.vcf --tmp_file ./tmp**
        cmdlist = ['sniffles', '-s', '3', '-z', '2', '-t', '4', '--report-seq', '--genotype','-m','merged.sorted.bam','-v','out.vcf']
        if subprocess.call(' '.join(cmdlist), shell=True)!= 0:
            print("sniffles sv fail")
            return
        print("sniffles sv is OK!")
    elif SoftwareSelector % 10 == 2:  # PBSV
        cmdlist = ['pbsv', 'discover', 'merged.sorted.bam', 'merged.sorted.svsig.gz']
        if subprocess.call(' '.join(cmdlist), shell=True)!= 0:
            print("pbsv discover is fail!")
            return
        print("pbsv discover is OK!")

        cmdlist = ['pbsv', 'call', Ref, 'merged.sorted.svsig.gz','out.vcf']
        if subprocess.call(' '.join(cmdlist), shell=True)!= 0:
            print("pbsv call is fail!")
            return
        print("pbsv call is OK!")

def RNA_snpindel(args):
    SeqRef=args.argR
    GTFRef=args.argG
    FastqFolder=args.argF
    SoftwareSelector=args.argS
    picardjaraddr=args.argP
    sampleName=args.argE
    print("RNA_NGS",SeqRef, GTFRef,FastqFolder,SoftwareSelector,picardjaraddr,sampleName)
    # SoftwareSelector
    # tens --- map 1:STAR    2:HISAT2
    # units --- variant 1:gatk

    SeqRefpath, _ = os.path.split(SeqRef)
    GTFRefpath, _ = os.path.split(GTFRef)
    # shotname, extension = os.path.splitext(tempfilename);
    CPUcount = cpu_count()
    mem = psutil.virtual_memory()
    TotalMem = int(mem.total) / 1024 / 1024 / 1024
    print(SeqRefpath,GTFRefpath, CPUcount, TotalMem)

    if SoftwareSelector // 10 == 1:  # STAR
        #index exist
        IndexExist = False
        if os.path.isfile('STARtmp/star_index/SAindex'):
            IndexExist = True

        # mkdir STARtmp dir
        cmdlist = ['mkdir', '-p', 'STARtmp/star_1pass', 'STARtmp/star_2pass', 'STARtmp/star_index','STARtmp/star_index_2pass']
        if subprocess.call(' '.join(cmdlist), shell=True)!= 0:
            print("mkdir is fail!")
            return
        print("mkdir is OK!")

        if not IndexExist:
            # STAR 1-pass index STAR --runThreadN 16 --runMode genomeGenerate --genomeDir ./STARtmp/star_index --genomeFastaFiles ../ref/hg19/human_g1k_v37.fasta --sjdbGTFfile ../ref/gtf/Homo_sapiens.GRCh37.75.gtf
            cmdlist = ['STAR', '--runThreadN', str(max(int(CPUcount*cpumemuse),4)), '--runMode', 'genomeGenerate','--genomeDir','STARtmp/star_index','--genomeFastaFiles',SeqRef,'--sjdbGTFfile',GTFRef]
            if subprocess.call(' '.join(cmdlist), shell=True)!= 0:
                print("1-pass index is fail!")
                return
            print("1-pass index is OK!")

        # STAR 1-pass align STAR --genomeDir
        cmdlist = ['STAR', '--runThreadN', str(max(int(CPUcount*cpumemuse),4)), '--genomeDir','STARtmp/star_index','--readFilesCommand','zcat','--outFileNamePrefix','STARtmp/star_1pass/my1pass','--readFilesIn']
        files = os.listdir(FastqFolder)
        for fqfile in files:
            if fqfile.endswith('gz'):
                cmdlist.extend([os.path.join(FastqFolder,fqfile)])
        print('STAR 1-pass: ', cmdlist)
        if subprocess.call(' '.join(cmdlist), shell=True)!= 0:
            print("1-pass align is fail!")
            return
        print("1-pass align is OK!")

        #STAR 2-pass index
        cmdlist = ['STAR', '--runThreadN', str(max(int(CPUcount * cpumemuse), 4)), '--runMode', 'genomeGenerate',
                   '--genomeDir', 'STARtmp/star_index_2pass', '--genomeFastaFiles', SeqRef, '--sjdbFileChrStartEnd ', "STARtmp/star_1pass/my1passSJ.out.tab"]
        if subprocess.call(' '.join(cmdlist), shell=True)!= 0:
            print("2-pass index is fail!")
            return
        print("2-pass index is OK!")

        #STAR 2-pass align
        cmdlist = ['STAR', '--runThreadN', str(max(int(CPUcount*cpumemuse),4)), '--genomeDir','STARtmp/star_index_2pass','--readFilesCommand','zcat','--outFileNamePrefix','STARtmp/star_2pass/my2pass','--readFilesIn']
        files = os.listdir(FastqFolder)
        for fqfile in files:
            if fqfile.endswith('gz'):
                cmdlist.extend([os.path.join(FastqFolder,fqfile)])
        print('STAR 2-pass align: ', cmdlist)
        if subprocess.call(' '.join(cmdlist), shell=True)!= 0:
            print("2-pass align is fail!")
            return
        print("2-pass align is OK!")

        # picard Add read groups, sort, mark duplicates, and create index
        cmdlist = ['java','-jar',picardjaraddr,'AddOrReplaceReadGroups','I=STARtmp/star_2pass/my2passAligned.out.sam','O=STARtmp/star_2pass/rg_added_sorted.bam','SO=coordinate','RGLB=rna','RGPL=illumina','RGPU=hiseq','RGSM='+sampleName]
        if subprocess.call(' '.join(cmdlist), shell=True)!= 0:
            print("picard AddOrReplaceReadGroups is fail!")
            return
        print("picard AddOrReplaceReadGroups is OK!")

        cmdlist = ['java','-jar',picardjaraddr,'MarkDuplicates','I=STARtmp/star_2pass/rg_added_sorted.bam','O={}.dedup.bam'.format(sampleName),'CREATE_INDEX=true','VALIDATION_STRINGENCY=SILENT','M=STARtmp/star_2pass/'+sampleName+'_dedup.metrics']
        if subprocess.call(' '.join(cmdlist), shell=True)!= 0:
            print("picard MarkDuplicates is fail!")
            return
        print("picard MarkDuplicates is OK!")
    elif SoftwareSelector // 10 == 2:  # HISAT2
        shotname, extension = os.path.splitext(GTFRef);
        # index exist
        IndexExist = True
        files = os.listdir(GTFRefpath)
        for fqfile in files:
            if fqfile.endswith('.ht2'):
                IndexExist=True
                break
        #create index
        if not IndexExist:
            cmdlist=['hisat2_extract_splice_sites.py',GTFRef,'>',os.path.join(GTFRefpath,shotname+'.ss') ]
            if subprocess.call(' '.join(cmdlist), shell=True) != 0:
                print("hisat2_extract_splice_sites is fail!")
                return
            print("hisat2_extract_splice_sites is OK!")
            cmdlist=['hisat2_extract_exons.py',GTFRef,'>',os.path.join(GTFRefpath,shotname+'.exon')]
            if subprocess.call(' '.join(cmdlist), shell=True) != 0:
                print("hisat2_extract_exons is fail!")
                return
            print("hisat2_extract_exons is OK!")
            cmdlist=['hisat2-build',SeqRef,'--exon',os.path.join(GTFRefpath,shotname+'.exon'),'-ss',os.path.join(GTFRefpath,shotname+'.ss'),'-p',str(max(int(CPUcount*cpumemuse),4)),sampleName]
            print("hisat2-build: ", cmdlist)
            if subprocess.call(' '.join(cmdlist), shell=True) != 0:
                print("hisat2-build is fail!")
                return
            print("hisat2-build is OK!")

        fqcount = 1
        files = os.listdir(FastqFolder)
        print(files)
        for fqfile in files:
            if fqfile.endswith('fq') or fqfile.endswith('fastq'):
                cmdlist = ['hisat2', '-x', sampleName, '-p', str(max(int(CPUcount*cpumemuse),4)), '-U',os.path.join(FastqFolder,fqfile), '-S',os.path.join(FastqFolder,'sub.{}.sam'.format(fqcount))]
                print('hisat2: ',cmdlist)
                if subprocess.call(' '.join(cmdlist), shell=True) != 0:
                    print("hisat2 is fail!")
                    return
                print("hisat2 is OK!")
                fqcount+=1

        # samtools sam---> sorted.bam
        files = os.listdir(FastqFolder)
        for fqfile in files:
            if fqfile.endswith('sam'):
                shotname, extension = os.path.splitext(fqfile)
                cmdlist = ['samtools', 'view', '-bS', os.path.join(FastqFolder,fqfile), '>', os.path.join(FastqFolder,shotname+'.bam')]
                if subprocess.call(' '.join(cmdlist), shell=True) != 0:
                    print("samtools: convert is fail")
                    return
        # merge all bam
        cmdlist=['samtools','merge','merged.sorted.bam']
        files = os.listdir(FastqFolder)
        for fqfile in files:
            if fqfile.endswith('bam'):
                cmdlist.extend([os.path.join(FastqFolder,fqfile)])
        if subprocess.call(' '.join(cmdlist), shell=True) != 0:
            print("samtools: merge is fail")
            return
        print("samtools: merge is OK!")


    if SoftwareSelector % 10 == 1: #gatk
        if SoftwareSelector //10 ==1:
            cmdlist=['gatk','SplitNCigarReads','-R',SeqRef,'-I','{}.dedup.bam'.format(sampleName),'-O','{}.dedup.split.bam'.format(sampleName)]
            if subprocess.call(' '.join(cmdlist), shell=True) != 0:
                print("gatk SplitNCigarReads is fail!")
                return
            print("gatk SplitNCigarReads is OK!")

        cmdlist = ['gatk', 'HaplotypeCaller ', '-R', SeqRef, '-I', '{}.dedup.split.bam'.format(sampleName), '-O', 'out.vcf']
        if subprocess.call(' '.join(cmdlist), shell=True) != 0:
            print("gatk HaplotypeCaller is fail!")
            return
        print("gatk HaplotypeCaller is OK!")

        # cmdlist = ['gatk', 'VariantFiltration ', '-R', SeqRef, '-V', 'out.vcf', '-window',
        #            '35', '-cluster', '3', "-filterName", "FS", "-filter", "FS > 30.0","-filterName", "QD", "-filter", "QD < 2.0" "-O",'out.filter.vcf']
        # if subprocess.call(' '.join(cmdlist), shell=True) != 0:
        #     print("gatk VariantFiltration is fail!")
        #     return
        # print("gatk VariantFiltration is OK!")

    print("ALL done!")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.description = 'DNAsv, RNAsnpindel'
    subparsers = parser.add_subparsers(help='sub-command help')

    # add DNA sv sub parser
    parser_DNAsv = subparsers.add_parser('DNAsv', help='DNAsv help')
    parser_DNAsv.add_argument("-R", help="human reference sequence", dest="argR", type=str)
    parser_DNAsv.add_argument("-F", help="fastq sequence folder", dest="argF", type=str)
    parser_DNAsv.add_argument("-S", help="software selector. DNA sv: tens --- map 1:minmap2  2:NGMLR   3:PBMM2, units --- variant 1:sniffles  2:PBSV. e.g. '21' means NGMLR+sniffles", dest="argS", type=str)
    # DNAsv default function
    parser_DNAsv.set_defaults(func=DNA_3GS)


    # add Rna snpindel sub parser
    parser_RNAsi = subparsers.add_parser('RNAsi', help='DNAsv help')
    parser_RNAsi.add_argument("-R", help="human reference sequence", dest="argR", type=str)
    parser_RNAsi.add_argument("-G", help="(RNA)human GTF reference, RNA-seq is required", dest="argG", type=str)
    parser_RNAsi.add_argument("-F", help="fastq sequence folder", dest="argF", type=str)
    parser_RNAsi.add_argument("-S", help="software selector. RNA snp&indel: tens---map 1-STAR 2-HISAT2; units---variant 1-gatk ", dest="argS", type=str)
    parser_RNAsi.add_argument("-P", help="picard.jar address", dest="argP", type=str, default='~/software/picard/picard.jar')
    parser_RNAsi.add_argument("-E", help="example name", dest="argE", type=str, default='mySample')
    # Rna snpindel default function
    parser_RNAsi.set_defaults(func=RNA_snpindel)

    args = parser.parse_args()
    args.func(args)
