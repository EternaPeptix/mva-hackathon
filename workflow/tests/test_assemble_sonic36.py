import copy
import tempfile
import unittest
from pathlib import Path

import assemble_sonic36 as m


class AssemblyControls(unittest.TestCase):
    def fixture(self):
        canonical=["BUB1B is unresolved.","A trial is proposed.","No result exists.","Review the evidence."]
        spoken=[m.spoken(s) for s in canonical]
        words=" ".join(spoken).split()
        return canonical,spoken,{"words":words,"start":[i*.2 for i in range(len(words))],"end":[i*.2+.15 for i in range(len(words))]}

    def test_canonical_preserved_and_spoken_expansion_aligned(self):
        c,s,t=self.fixture();rows=m.align_sentences(c,s,t,10)
        self.assertEqual([r['caption'] for r in rows],c)
        self.assertEqual(rows[0]['word_end_index_exclusive'],5)

    def test_word_tamper_refused(self):
        c,s,t=self.fixture();t['words'][1]='two'
        with self.assertRaisesRegex(ValueError,'word stream'):m.align_sentences(c,s,t,10)

    def test_missing_word_refused(self):
        c,s,t=self.fixture();t['words'].pop()
        with self.assertRaises(ValueError):m.align_sentences(c,s,t,10)

    def test_changed_canonical_expansion_refused(self):
        c,s,t=self.fixture();s[0]='A new claim.'
        with self.assertRaises(ValueError):m.align_sentences(c,s,t,10)

    def test_timestamp_bounds_refused(self):
        c,s,t=self.fixture();t['end'][-1]=11
        with self.assertRaisesRegex(ValueError,'outside'):m.align_sentences(c,s,t,10)

    def test_nonfinite_timestamp_refused(self):
        c,s,t=self.fixture();t['start'][0]=float('nan')
        with self.assertRaises(ValueError):m.align_sentences(c,s,t,10)

    def test_backward_timing_refused(self):
        c,s,t=self.fixture();t['start'][3]=0
        with self.assertRaisesRegex(ValueError,'Nonmonotone'):m.align_sentences(c,s,t,10)

    def test_zero_duration_word_retained(self):
        c,s,t=self.fixture();t['end'][1]=t['start'][1]
        self.assertEqual(len(m.align_sentences(c,s,t,10)),4)

    def test_provider_segment_inside_sentence_supported(self):
        c,s,t=self.fixture();t['words'][0]=' '.join(t['words'][:2]);t['words'].pop(1);t['start'].pop(1);t['end'].pop(0)
        self.assertEqual(len(m.align_sentences(c,s,t,10)),4)

    def test_provider_segment_crosses_caption_refused(self):
        c,s,t=self.fixture();t['words'][4]+=' '+t['words'][5];t['words'].pop(5);t['start'].pop(5);t['end'].pop(4)
        with self.assertRaisesRegex(ValueError,'crosses'):m.align_sentences(c,s,t,10)

    def test_continuous_scene_lengths_and_final_boundaries(self):
        counts=[round(n*m.SR) for n in [18.56,20.24,22.88,23.28,27.2,33.12]]
        boundaries,pad=m.frame_boundaries(counts)
        self.assertEqual((boundaries[0],boundaries[-1]),(0,4320));self.assertAlmostEqual(pad,34.72)
        self.assertEqual(len(boundaries),7)

    def test_overlong_audio_refused_not_trimmed(self):
        with self.assertRaises(ValueError):m.frame_boundaries([31*m.SR]*6)

    def test_excessive_reading_pad_refused(self):
        with self.assertRaises(ValueError):m.frame_boundaries([10*m.SR]*6)

    def test_stale_source_pin_refused(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'source';p.write_text('original');info=m.pin(p);p.write_text('modified')
            with self.assertRaises(ValueError):m.bind(info)


if __name__=='__main__':unittest.main()
