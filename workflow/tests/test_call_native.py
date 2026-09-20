"""Real native callers must recover held-out synthetic SNV/indel truth and resume safely."""
import gzip
import importlib.util
import json
import os
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest

import pysam

spec = importlib.util.spec_from_file_location('call_genome', Path(__file__).parents[1]/'scripts/call_genome.py')
caller = importlib.util.module_from_spec(spec)
spec.loader.exec_module(caller)


class ShardTests(unittest.TestCase):
    def test_complete_contigs_once_in_reference_order(self):
        contigs=[('1',300),('2',200),('X',10),('decoy1',2),('decoy2',3)]
        groups=caller.shard_contigs(contigs,100,2)
        self.assertEqual([x for group in groups for x in group],contigs)
        self.assertEqual(groups,[contigs[:1],contigs[1:2],contigs[2:4],contigs[4:]])
        with self.assertRaises(ValueError):
            caller.shard_contigs([('broken',0)])


class NativeCallerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.root=Path(os.environ.get('MVA_PROJECT_ROOT','/Users/jeweled/Documents/mva-hackathon'))
        cls.bin=cls.root/'tools/bio/bin'
        cls.jar=cls.root/'gatk-4.6.1.0/gatk-package-4.6.1.0-local.jar'
        if not cls.jar.exists() or not (cls.bin/'bwa').exists():
            raise unittest.SkipTest('native Studio runtime required')

    def test_two_callers_recover_snv_and_indel_and_preserve_completed_run(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder=Path(tmp)
            rng=random.Random(20260919)
            sequences={c:''.join(rng.choice('ACGT') for _ in range(3000)) for c in ('1','2')}
            snv=1000
            alt='A' if sequences['1'][snv]!='A' else 'C'
            deletion=1000
            while sequences['2'][deletion]==sequences['2'][deletion+2]: deletion+=1
            mutated={'1':sequences['1'][:snv]+alt+sequences['1'][snv+1:],
                     '2':sequences['2'][:deletion+1]+sequences['2'][deletion+3:]}
            truth={('1',snv+1,sequences['1'][snv],alt),
                   ('2',deletion+1,sequences['2'][deletion:deletion+3],sequences['2'][deletion])}
            reference=folder/'reference.fa'
            reference.write_text(''.join('>'+c+'\n'+s+'\n' for c,s in sequences.items()))
            pysam.faidx(str(reference))
            subprocess.run([str(self.bin/'bwa'),'index',str(reference)],capture_output=True,check=True)
            files=[folder/'R1.fastq.gz',folder/'R2.fastq.gz']
            for mate,p in enumerate(files,1):
                with gzip.open(p,'wt') as fastq:
                    for contig in sequences:
                        for allele in (0,1):
                            genome=(sequences if allele==0 else mutated)[contig]
                            for index in range(24):
                                start=900+index*2+(220 if mate==2 else 0)
                                read=genome[start:start+150]
                                if mate==2: read=read.translate(str.maketrans('ACGT','TGCA'))[::-1]
                                name=f'SIM:{contig}:{allele}:{index}'
                                fastq.write('@'+name+f' {mate}:N:0:ACGT\n'+read+'\n+\n'+'I'*len(read)+'\n')
            sam=folder/'reads.sam'
            with sam.open('w') as stream:
                subprocess.run([str(self.bin/'bwa'),'mem','-t','2','-R',r'@RG\tID:SIM\tSM:SYNTH\tPL:ILLUMINA',str(reference),*[str(p) for p in files]],stdout=stream,stderr=subprocess.PIPE,check=True)
            bam=folder/'reads.bam'
            subprocess.run([str(self.bin/'samtools'),'sort','-o',str(bam),str(sam)],capture_output=True,check=True)
            subprocess.run([str(self.bin/'samtools'),'index',str(bam)],capture_output=True,check=True)
            output=folder/'calls'
            command=[sys.executable,str(Path(caller.__file__)),'--root',str(self.root),'--reference',str(reference),
                     '--alignment',str(bam),'--sample','SYNTH','--output',str(output),'--workers','2','--shard-max-bases','4000']
            result=subprocess.run(command,capture_output=True,text=True)
            if result.returncode:
                logs='\n'.join(p.read_text(errors='replace')[-2000:] for p in output.rglob('*.log'))
                self.fail(result.stdout+result.stderr+logs)
            status=json.loads((output/'status.json').read_text())
            self.assertEqual(status['status'],'complete')
            self.assertEqual(status['expected_shards'],2)
            mtimes={}
            for method,data in status['outputs'].items():
                path=Path(data['path']);mtimes[method]=path.stat().st_mtime_ns
                with pysam.VariantFile(str(path)) as vcf:
                    found={(r.chrom,r.pos,r.ref,r.alts[0]):r.samples['SYNTH']['GT'] for r in vcf.fetch()}
                for site in truth:
                    self.assertIn(site,found,(method,found,truth))
                    self.assertEqual(tuple(sorted(found[site])),(0,1),(method,site,found[site]))
            repeat=subprocess.run(command,capture_output=True,text=True)
            self.assertEqual(repeat.returncode,0,repeat.stdout+repeat.stderr)
            for method,data in status['outputs'].items():
                self.assertEqual(Path(data['path']).stat().st_mtime_ns,mtimes[method])
            # A different method/input definition must never silently reuse old calls.
            changed=subprocess.run(command[:-1]+['9000'],capture_output=True,text=True)
            self.assertNotEqual(changed.returncode,0)


if __name__=='__main__':
    unittest.main()
