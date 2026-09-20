#!/usr/bin/env python3
"""Create original schematic slides and sentence-aligned local narration/captions.

All visuals are explanatory diagrams, never simulated patient measurements.
Generation runs on the Mac Studio using an installed system voice and Pillow.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import textwrap
import time
import wave

import numpy as np
from PIL import Image, ImageDraw, ImageFont

W,H,FPS,SR=1920,1080,24,44100
BG='#F4F6F2';INK='#142F37';MUTED='#4E676D';TEAL='#087F80';PURPLE='#6D57A0';AMBER='#A66620';PALE='#E3ECE8';WHITE='#FFFFFF'
FONT=Path('/System/Library/Fonts/Supplemental/Arial.ttf')
BOLD=Path('/System/Library/Fonts/Supplemental/Arial Bold.ttf')


def save(path,data):path.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def font(size,bold=False):return ImageFont.truetype(str(BOLD if bold else FONT),size)


class Slide:
    def __init__(self,scene,index,footer):
        self.image=Image.new('RGB',(W,H),BG);self.d=ImageDraw.Draw(self.image);self.bounds=[]
        self.text((96,54),'JEWELED  /  MVA RESEARCH',font(24,True),MUTED)
        self.d.rounded_rectangle((1458,38,1824,88),radius=24,fill=PALE)
        self.text((1484,52),'DRAFT   ·   '+str(index+1)+' / 6',font(23,True),TEAL)
        self.text((96,132),scene['kicker'],font(25,True),TEAL)
        self.text((96,182),scene['title'],font(62,True),INK,max_width=1728)
        self.text((96,808),scene['display_source'],font(22),MUTED,max_width=1728)
        self.text((96,1035),footer,font(20),MUTED,max_width=1728)

    def text(self,xy,value,f,fill=INK,max_width=None,line_spacing=8):
        lines=[]
        for paragraph in value.split('\n'):
            if not max_width:lines.append(paragraph);continue
            line=''
            for word in paragraph.split():
                candidate=(line+' '+word).strip()
                if self.d.textlength(candidate,font=f)>max_width:
                    if not line:raise ValueError('word does not fit')
                    lines.append(line);line=word
                else:line=candidate
            lines.append(line)
        text='\n'.join(lines)
        bbox=self.d.multiline_textbbox(xy,text,font=f,spacing=line_spacing)
        if min(bbox[:2])<0 or bbox[2]>W or bbox[3]>H:raise ValueError('text outside frame: '+value)
        self.d.multiline_text(xy,text,font=f,fill=fill,spacing=line_spacing)
        self.bounds.append({'text':value,'box':list(bbox)})
        return bbox

    def card(self,box,title,body,color=TEAL):
        self.d.rounded_rectangle(box,radius=24,fill=WHITE,outline='#CEDBD6',width=2)
        x,y,x2,y2=box
        self.d.rounded_rectangle((x+28,y+30,x+38,y2-30),radius=4,fill=color)
        b=self.text((x+66,y+40),title,font(42,True),color,max_width=x2-x-102)
        b=self.text((x+66,b[3]+28),body,font(34),INK,max_width=x2-x-108,line_spacing=13)
        if b[3]>y2-25:raise ValueError('card text overflow')

    def arrow(self,x1,y1,x2,y2,color=TEAL,width=6):
        self.d.line((x1,y1,x2,y2),fill=color,width=width)
        angle=math.atan2(y2-y1,x2-x1);length=19
        points=[(x2,y2),(x2-length*math.cos(angle-.55),y2-length*math.sin(angle-.55)),(x2-length*math.cos(angle+.55),y2-length*math.sin(angle+.55))]
        self.d.polygon(points,fill=color)


def make_slide(scene,index,footer):
    s=Slide(scene,index,footer);kind=scene['id']
    if kind=='purpose':
        s.card((96,318,636,724),'Genome','Which explanation fits the evidence?')
        s.card((692,318,1232,724),'Cell function','Which defect can be measured and corrected?')
        s.card((1288,318,1824,724),'Next decision','Advance, revise or end a hypothesis.')
        s.arrow(648,518,680,518);s.arrow(1244,518,1276,518)
    elif kind=='genotype':
        for offset,label,opposite in [(96,'Possible: different copies',True),(990,'Possible: same copy',False)]:
            s.d.rounded_rectangle((offset,310,offset+834,740),radius=24,fill=WHITE,outline='#CEDBD6',width=2)
            s.text((offset+34,343),label,font(38,True),TEAL)
            for n,y in enumerate((480,640)):
                s.text((offset+35,y-57),'Copy '+str(n+1),font(25),MUTED)
                s.d.rounded_rectangle((offset+36,y-10,offset+788,y+10),radius=10,fill='#C8D9D3')
            marks=[(offset+250,480,'Stop',AMBER),(offset+586,640 if opposite else 480,'Missense',PURPLE)]
            for x,y,label,color in marks:
                s.d.ellipse((x-17,y-17,x+17,y+17),fill=color)
                width=s.d.textlength(label,font=font(30,True))
                s.text((x-width/2,y+30),label,font(30,True),color)
        s.text((96,759),'Phase remains unknown. Positions are schematic.',font(27,True),MUTED)
    elif kind=='mechanism':
        s.card((96,315,924,700),'Restore mitotic function','Protein quantity / quality\n↓\nAccurate chromosome segregation',TEAL)
        s.card((996,315,1824,700),'Protect downstream function','Demonstrated cellular stress\n↓\nFunction and viable-cell output',PURPLE)
        s.text((96,746),'Start with a reproducible defect and an interpretable genetic rescue control.',font(30,True),INK,max_width=1728)
    elif kind=='candidates':
        s.card((96,315,924,680),'Arimoclomol','Conditional protein-quality probe\n\nBUBR1 rescue is unproven.',TEAL)
        s.card((996,315,1824,680),'Acetylcysteine','Exploratory stress comparator\n\nGenotype-specific benefit is unproven.',PURPLE)
        s.text((96,732),'Test separately  →  Measure exposure  →  Require functional improvement',font(34,True),INK,max_width=1728)
    elif kind=='measurement':
        s.card((96,310,636,647),'Division','Errors ÷ observed divisions',TEAL)
        s.card((692,310,1232,647),'Function','An independent functional readout',PURPLE)
        s.card((1288,310,1824,647),'Cell output','Starting cells, viable cells and lineage',TEAL)
        s.d.rounded_rectangle((96,690,1824,774),radius=18,fill='#F0E5D7')
        s.text((130,711),'Check for arrest, selective death and abnormal-clone expansion.',font(32,True),AMBER,max_width=1650)
    elif kind=='decision':
        rows=[('01','Finish genomic review and establish model feasibility.'),('02','Demonstrate function at justified exposure.'),('03','Confirm independently and review the translational case.')]
        for i,(number,body) in enumerate(rows):
            y=314+i*142;s.d.rounded_rectangle((96,y,1824,y+112),radius=20,fill=WHITE,outline='#CEDBD6',width=2)
            s.text((132,y+31),number,font(42,True),TEAL)
            s.text((248,y+35),body,font(34,True),INK,max_width=1510)
        s.text((96,763),'Laboratory access and drug response are not established.',font(27,True),MUTED)
    else:raise ValueError('unknown scene')
    return s


def run_command(args,timeout=90):
    result=subprocess.run(list(map(str,args)),capture_output=True,text=True,timeout=timeout)
    if result.returncode:raise RuntimeError(str(args[0])+' failed: '+result.stderr)
    return {'command':list(map(str,args)),'exit_code':result.returncode}


def synthesize(story,folder,rate):
    folder.mkdir();clips=[];commands=[]
    for i,scene in enumerate(story['scenes']):
        for j,sentence in enumerate(scene['sentences']):
            stem=f'{i+1:02d}-{j+1:02d}';text=folder/(stem+'.txt');aiff=folder/(stem+'.aiff');wav=folder/(stem+'.wav')
            # Expand gene/protein names for the generic system voice; captions retain canonical names.
            spoken=sentence.replace('BUB1B','Bub one B').replace('BUBR1','Bub R one')
            text.write_text(spoken+'\n')
            commands.append(run_command(['/usr/bin/say','-v',story['synthetic_narration']['voice'],'-r',str(rate),'-f',text,'-o',aiff]))
            commands.append(run_command(['/usr/bin/afconvert',aiff,wav,'-f','WAVE','-d','LEI16@44100','-c','1']))
            with wave.open(str(wav),'rb') as f:
                if (f.getframerate(),f.getnchannels(),f.getsampwidth())!=(SR,1,2):raise ValueError('unexpected speech format')
                count=f.getnframes();data=np.frombuffer(f.readframes(count),dtype='<i2')
            if not len(data) or not np.any(data):raise ValueError('empty/silent speech clip')
            clips.append({'scene':i,'sentence':j,'caption':sentence,'spoken_text':spoken,'wav':str(wav),'samples':len(data),'seconds':len(data)/SR,'sha256':sha(wav)})
    save(folder/'commands.json',commands);save(folder/'clips.json',clips);return clips


def caption_time(t):
    ms=round(t*1000);hours,ms=divmod(ms,3600000);minutes,ms=divmod(ms,60000);seconds,ms=divmod(ms,1000)
    return f'{hours:02d}:{minutes:02d}:{seconds:02d}.{ms:03d}'


def run(story_path,output):
    output.mkdir(parents=True,exist_ok=False);tick=time.monotonic();os.nice(10)
    story=json.loads(story_path.read_text());assert story['duration_seconds']==180 and len(story['scenes'])==6
    rate=story['synthetic_narration']['initial_words_per_minute'];attempts=[]
    while True:
        clips=synthesize(story,output/f'voice-{rate}',rate)
        duration=sum(c['seconds'] for c in clips)+.18*(len(clips)-6)+.7*6
        attempts.append({'rate':rate,'speech_plus_minimum_pauses_seconds':duration})
        if duration<=174:break
        rate=math.ceil(rate*duration/171)
        if rate>170:raise ValueError('narration too long at a comfortable rate; revise text instead of truncating')
    scenes=[];base_images=[];bounds=[]
    for i,scene in enumerate(story['scenes']):
        slide=make_slide(scene,i,story['persistent_footer']);p=output/f'scene-{i+1:02d}.png';slide.image.save(p)
        base_images.append(slide.image);bounds.extend(slide.bounds);scenes.append({'id':scene['id'],'image':str(p),'title':scene['title']})
    # Give each slide equal extra reading time after accounting for its actual speech.
    per_scene=[sum(c['seconds'] for c in clips if c['scene']==i)+.18*3+.7 for i in range(6)]
    extra=(180-sum(per_scene))/6
    if extra<1:raise ValueError('insufficient reading time')
    boundaries=[0]
    for seconds in per_scene:boundaries.append(boundaries[-1]+(seconds+extra)*FPS)
    frame_boundaries=[round(x) for x in boundaries];frame_boundaries[-1]=180*FPS
    narration=np.zeros(180*SR,dtype=np.float64);events=[];captions=[];segments=[]
    for i,scene in enumerate(story['scenes']):
        start_frame,end_frame=frame_boundaries[i:i+2];cursor=start_frame/FPS+.7;last_frame=start_frame
        for c in [x for x in clips if x['scene']==i]:
            begin=round(cursor*SR);end=begin+c['samples']
            if end>round(end_frame/FPS*SR):raise ValueError('speech exceeds its scene')
            with wave.open(c['wav'],'rb') as f:data=np.frombuffer(f.readframes(c['samples']),dtype='<i2').astype(float)
            if np.any(narration[begin:end]):raise ValueError('speech overlaps')
            narration[begin:end]=data
            cstart=round(cursor*FPS);cend=math.ceil((cursor+c['seconds'])*FPS)
            if cstart>last_frame:segments.append({'start_frame':last_frame,'end_frame':cstart,'image':scenes[i]['image']})
            image=base_images[i].copy();d=ImageDraw.Draw(image);f=font(38);lines=[];line=''
            for word in c['caption'].split():
                candidate=(line+' '+word).strip()
                if d.textlength(candidate,font=f)>1636:lines.append(line);line=word
                else:line=candidate
            lines.append(line)
            if len(lines)>3:raise ValueError('caption too dense')
            d.rounded_rectangle((96,856,1824,1005),radius=18,fill=INK)
            text='\n'.join(lines);box=d.multiline_textbbox((140,874),text,font=f,spacing=9)
            if box[2]>1784 or box[3]>995:raise ValueError('caption overflow')
            d.multiline_text((140,874),text,font=f,fill=WHITE,spacing=9)
            p=output/f'caption-{len(captions)+1:02d}.png';image.save(p)
            segments.append({'start_frame':cstart,'end_frame':cend,'image':str(p)})
            event={**c,'start_seconds':begin/SR,'end_seconds':end/SR,'start_sample':begin,'end_sample':end,'caption_image':str(p)}
            captions.append(event);last_frame=cend;cursor=end/SR+.18
        if last_frame<end_frame:segments.append({'start_frame':last_frame,'end_frame':end_frame,'image':scenes[i]['image']})
        scenes[i].update(start_frame=start_frame,end_frame=end_frame,start_seconds=start_frame/FPS,end_seconds=end_frame/FPS)
    if segments[0]['start_frame']!=0 or segments[-1]['end_frame']!=4320:raise ValueError('video endpoints incomplete')
    for a,b in zip(segments,segments[1:]):
        if a['end_frame']!=b['start_frame']:raise ValueError('video frame gap/overlap')
    peak=float(np.max(np.abs(narration)));gain=(32767*10**(-3/20))/peak
    encoded=np.round(narration*gain).astype('<i2');audio=output/'narration.wav'
    with wave.open(str(audio),'wb') as w:w.setnchannels(1);w.setsampwidth(2);w.setframerate(SR);w.writeframes(encoded.tobytes())
    vtt='WEBVTT\n\n'+''.join(f"{n+1}\n{caption_time(c['start_seconds'])} --> {caption_time(c['end_seconds'])}\n{c['caption']}\n\n" for n,c in enumerate(captions))
    (output/'captions.vtt').write_text(vtt)
    sheet=Image.new('RGB',(1280,1080),WHITE)
    for i,img in enumerate(base_images):sheet.paste(img.resize((640,360)),((i%2)*640,(i//2)*360))
    sheet.save(output/'contact-sheet.png')
    plan={'title':story['title'],'draft':True,'width':W,'height':H,'fps':FPS,'duration_seconds':180,'frames':4320,'audio':str(audio),'segments':segments,'scenes':scenes,'captions':captions}
    save(output/'video-plan.json',plan)
    receipt={'status':'complete','elapsed_seconds':time.monotonic()-tick,'storyboard_sha256':sha(story_path),'code_sha256':sha(Path(__file__)),'voice':story['synthetic_narration']['voice'],'voice_attempts':attempts,'selected_rate':rate,'captions':len(captions),'narrated_words':sum(len(c['caption'].split()) for c in captions),'sample_rate':SR,'audio_samples':len(encoded),'audio_seconds':len(encoded)/SR,'peak_dBFS':20*math.log10(float(np.max(np.abs(encoded.astype(float))))/32767),'global_audio_gain':gain,'fonts':{str(p):sha(p) for p in [FONT,BOLD]},'diagram_text_bounds':bounds,'files':{p.name:{'bytes':p.stat().st_size,'sha256':sha(p)} for p in output.iterdir() if p.is_file()},'limits':['Research draft; complete genome analysis and laboratory feasibility pending','Schematic illustrations contain no experimental data','Generic system voice, not the user or a patient','A local render is not a published submission']}
    save(output/'preparation.json',receipt);return receipt


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ('storyboard','output'):parser.add_argument('--'+key,required=True,type=Path)
    args=parser.parse_args();r=run(args.storyboard,args.output);print(json.dumps({k:r[k] for k in ['status','elapsed_seconds','selected_rate','captions','narrated_words','audio_seconds','peak_dBFS']}))
