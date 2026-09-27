"""Empirical injections: fixed populations, grouped photons and label controls."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

HERE = Path(__file__).resolve().parents[1] / 'studies/synthetic-timing-shifts'


def module(name, filename):
    spec = importlib.util.spec_from_file_location(name, HERE / filename)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


empirical = module('empirical_timing_test', 'empirical.py')
extractor = module('empirical_extract_test', 'empirical_extract.py')


def background():
    # Several local timing shapes, four strata, two chronological segments.
    group = np.repeat(np.arange(8), 40)
    return dict(times=np.tile(np.linspace(-3, 4, 40), 8), fold=group//4,
                strata=group % 4, regime=(group % 4)//2,
                population=np.full((2, 4), 10000, dtype=np.int64),
                event_counts=np.full((2, 4), 40, dtype=np.int64))


class EmpiricalTimingTests(unittest.TestCase):
    def test_artificial_assignments_preserve_all_trial_exposures(self):
        b = background()
        x, exposures = empirical.assign(b, np.random.default_rng(7))
        np.testing.assert_array_equal(exposures.sum(axis=2), b['population'])
        self.assertEqual(len(x), int(b['event_counts'].sum()))
        self.assertTrue(np.all(exposures >= 0))

    def test_every_shift_retains_event_population_and_count_score(self):
        b = background()
        x, exposures = empirical.assign(b, np.random.default_rng(13))
        baseline = empirical.histogram(b, x, 0, 'fixed')
        for mode, shift in [('fixed', .125), ('reversing', 100),
                            ('transfer_failure', 1), ('local_only', 1),
                            ('digital_fixed', 2), ('digital_reversing', 4),
                            ('redigitized_fixed', .25), ('redigitized_reversing', .5)]:
            got = empirical.histogram(b, x, shift, mode, np.zeros(len(x)))
            np.testing.assert_array_equal(got.sum(axis=-1), baseline.sum(axis=-1))
            self.assertEqual(got.sum(), len(x))
            before = empirical.test_scores(baseline[1], empirical.learn(baseline[0], exposures[0]))
            after = empirical.test_scores(got[1], empirical.learn(got[0], exposures[0]))
            self.assertEqual(before['click_count'], after['click_count'])

    def test_zero_shift_and_zero_redigitized_shift_are_identical(self):
        b = background()
        rng = np.random.default_rng(23)
        x, _ = empirical.assign(b, rng)
        latent = rng.uniform(-.5, .5, len(x))
        baseline = empirical.histogram(b, x, 0, 'fixed')
        for mode in ('reversing', 'local_only', 'transfer_failure', 'digital_fixed',
                     'redigitized_fixed', 'redigitized_reversing'):
            np.testing.assert_array_equal(baseline, empirical.histogram(b, x, 0, mode, latent))

    def test_redigitization_counts_integer_crossings(self):
        # A +0.25-bin arm shift crosses exactly one quarter of uniform latent
        # positions; the opposite arm crosses in the opposite direction.
        b = dict(times=np.full(8, .25), fold=np.zeros(8, dtype=int),
                 regime=np.zeros(8, dtype=int), strata=np.zeros(8, dtype=int))
        latent = np.tile(np.array([-.375, -.125, .125, .375]), 2)
        x = np.repeat([0, 1], 4)
        got = empirical.histogram(b, x, .5, 'redigitized_fixed', latent)
        self.assertEqual(got[0, 0, 0, 8], 1)  # one negative crossing
        self.assertEqual(got[0, 0, 0, 9], 3)
        self.assertEqual(got[0, 0, 1, 9], 3)
        self.assertEqual(got[0, 0, 1, 10], 1)  # one positive crossing

    def test_learning_uses_separate_training_counts(self):
        counts = np.zeros((4, 2, 18), dtype=np.int64)
        for state in range(4):
            counts[state, 0, 7] = 50
            counts[state, 1, 10] = 50
        exposures = np.full((4, 2), 1000)
        features = empirical.learn(counts, exposures)
        before = {k: v.copy() for k, v in features.items()}
        empirical.test_scores(counts[:, ::-1].copy(), features)
        for name in features:
            np.testing.assert_array_equal(features[name], before[name])
        self.assertLess(empirical.test_scores(counts, features)['trained_timing'], 1e-10)

    def test_unobserved_centered_latent_model_can_hide_small_physical_shifts(self):
        b = background()
        x, _ = empirical.assign(b, np.random.default_rng(54))
        null = empirical.histogram(b, x, 0, 'fixed')
        hidden = empirical.histogram(b, x, .25, 'centered_latent_reversing')
        np.testing.assert_array_equal(hidden, null)

    def test_multiclicks_are_retained_but_receive_one_trial_label(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)/'cache.npz'
            population = np.full((32, 2), 5, dtype=np.int64)
            provenance = dict(selected_events=4, retained_trials=int(population.sum()))
            np.savez(path, event_row=np.array([0, 0, 10, 10]),
                     event_phase=np.array([.1, 3., .2, -4.]), event_state=np.array([0, 0, 1, 1]),
                     source_rows=np.array(320), exposures=population,
                     provenance=np.array(json.dumps(provenance)))
            b = empirical.load_background(path)
            np.testing.assert_array_equal(b['times'], [.1, .2])
            self.assertEqual(b['event_counts'].sum(), 2)
            x, exposures = empirical.assign(b, np.random.default_rng(4))
            self.assertEqual(len(x), 2)
            self.assertEqual(exposures.sum(), 320)

    def test_optional_raw_timing_does_not_change_existing_decoder(self):
        decoder = extractor.decoder
        config = dict(pk=125, radius=5, bitoffset=37)
        a = np.array([(6, 0, 0), (2, 1, 0), (0, 6741, 0), (0, 6902, 0),
                      (6, 129102, 0), (4, 129103, 0), (0, 135843, 0)], dtype=decoder.DTYPE)
        old = decoder.decode_block(a, 258204, None, config)
        new = decoder.decode_block(a, 258204, None, config, include_events=True)
        np.testing.assert_array_equal(old[0], new[0])
        np.testing.assert_array_equal(old[2], new[2])
        self.assertEqual(old[3:], new[3:5])
        for name in decoder.VARIANTS:
            np.testing.assert_array_equal(old[1][name], new[1][name])
        events = new[-1]
        np.testing.assert_array_equal(events['row'], [0, 0, 1])
        recovered = np.zeros(2, dtype='u2')
        bit = events['pulse']-config['bitoffset']
        eligible = (np.abs(events['phase']-config['pk']) < config['radius']) & (bit >= 0) & (bit < 16)
        np.bitwise_or.at(recovered, events['row'][eligible], (1 << bit[eligible]).astype('u2'))
        np.testing.assert_array_equal(recovered, old[1]['nominal'])

    def test_refuse_in_repo_event_cache_and_invalid_injection(self):
        with self.assertRaises(ValueError):
            extractor.outside_repository(HERE/'source.npz')
        b = background()
        x = np.zeros(len(b['times']), dtype=int)
        for shift, mode in [(1, 'digital_fixed'), (1, 'unknown'), (.1, 'redigitized_fixed')]:
            with self.assertRaises(ValueError):
                empirical.histogram(b, x, shift, mode)


if __name__ == '__main__':
    unittest.main()
