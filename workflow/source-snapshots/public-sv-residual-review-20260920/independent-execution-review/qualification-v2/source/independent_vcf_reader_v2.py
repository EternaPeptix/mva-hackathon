"""Read saved VCF text only; no producer imports, index, alignment or subprocess.

This helper preserves literal strings and INFO presence independently of typed
values. Unselected original-context rows need not be retained in memory.
"""
import gzip
import math
from pathlib import Path

ANNOTATIONS=('TRF','TRFdiff','TRFrepeat','TRFovl','TRFstart','TRFend','TRFperiod','TRFcopies','TRFscore','TRFentropy','TRFsim','RM_clsfam','RM_repeat','RM_score','LCR','REMAP')

def declaration(line):
    """Split a structured header outside quoted/escaped text."""
    text=line[line.index('<')+1:-1]
    cells=[];buf=[];quoted=False;escaped=False
    for char in text:
        if escaped:buf.append(char);escaped=False;continue
        if char=='\\' and quoted:buf.append(char);escaped=True;continue
        if char=='"':quoted=not quoted;buf.append(char);continue
        if char==',' and not quoted:cells.append(''.join(buf));buf=[]
        else:buf.append(char)
    if quoted or escaped:raise ValueError('Unterminated structured header')
    cells.append(''.join(buf));out={}
    for cell in cells:
        name,sep,value=cell.partition('=')
        if not sep or name in out:raise ValueError('Malformed/duplicate header field')
        out[name]=value
    if not all(k in out for k in ('ID','Number','Type')):raise ValueError('Incomplete INFO declaration')
    return out


def info_fields(text):
    if text=='.':return {}
    out={}
    for cell in text.split(';'):
        name,sep,value=cell.partition('=')
        if not name or name in out:raise ValueError('Empty/duplicate INFO key')
        out[name]=value if sep else None
    return out


def annotation(info,tag,definitions):
    if tag not in definitions:raise ValueError('Selected INFO header absent: '+tag)
    decl=definitions[tag];kind=decl['Type'];number=decl['Number'];present=tag in info
    raw=info.get(tag)
    result={'header_number':number,'header_type':kind,'present':present,'raw':raw,'value':None}
    if not present:return result
    if kind=='Flag':
        if number!='0' or raw is not None:raise ValueError('Flag must be valueless Number=0')
        result['value']=True;return result
    if raw is None:raise ValueError('Nonflag INFO lacks value')
    parts=raw.split(',')
    if number.isdecimal() and len(parts)!=int(number):raise ValueError('INFO cardinality disagrees')
    def scalar(s):
        if s=='.':return None
        if kind=='Integer':
            value=int(s)
            if not s or any(c not in '+-0123456789' for c in s):raise ValueError('Invalid integer')
            return value
        if kind=='Float':
            value=float(s)
            if not math.isfinite(value):raise ValueError('Nonfinite numeric INFO')
            return value
        if kind in ('String','Character'):
            if kind=='Character' and len(s)!=1:raise ValueError('Invalid Character INFO')
            return s
        raise ValueError('Unsupported INFO type')
    vals=[scalar(s) for s in parts]
    result['value']=vals[0] if number=='1' else vals
    return result


def identity(row):return (row['chrom'],row['pos'],row['id'],row['ref'],tuple(row['alts']))


def read_vcf(path,keep=None):
    """Return header, retained rows and full count; keep is an identity set or None."""
    p=Path(path);opener=gzip.open if p.suffix=='.gz' else open
    headers=[];defs={};samples=None;rows=[];count=0
    with opener(p,'rt',encoding='utf-8',newline='') as handle:
        for line in handle:
            text=line.rstrip('\n\r')
            if text.startswith('##'):
                if samples is not None:raise ValueError('Metadata after column header')
                headers.append(text)
                if text.startswith('##INFO=<'):
                    d=declaration(text)
                    if d['ID'] in defs:raise ValueError('Duplicate INFO declaration')
                    defs[d['ID']]=d
                continue
            if text.startswith('#CHROM\t'):
                if samples is not None:raise ValueError('Duplicate column header')
                fields=text.split('\t')
                if fields[:9]!=['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']:raise ValueError('Unsupported fixed header')
                samples=fields[9:]
                if len(samples)!=1:raise ValueError('Require exactly one sample')
                headers.append(text);continue
            if samples is None or not text:raise ValueError('VCF data/header ordering')
            cells=text.split('\t')
            if len(cells)!=10:raise ValueError('VCF row/sample count')
            count+=1
            row={'ordinal':count,'chrom':cells[0],'pos':int(cells[1]),'id':None if cells[2]=='.' else cells[2],
                 'ref':cells[3],'alts':cells[4].split(','),'qual_literal':cells[5],'filter_literal':cells[6],
                 'info_literal':cells[7],'format_literal':cells[8],'sample_literal':cells[9],'record_literal':text}
            if keep is not None and identity(row) not in keep:continue
            fmt=cells[8].split(':');vals=cells[9].split(':')
            if len(set(fmt))!=len(fmt) or len(vals)>len(fmt):raise ValueError('FORMAT duplicates or value count')
            row['format_values']=dict(zip(fmt,vals));row['info']=info_fields(cells[7])
            row['annotations']={tag:annotation(row['info'],tag,defs) for tag in ANNOTATIONS}
            rows.append(row)
    if samples is None:raise ValueError('No column header')
    return {'headers':headers,'info_definitions':defs,'samples':samples,'total_records':count,'rows':rows}
