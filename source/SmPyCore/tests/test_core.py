"""Focused smoke tests for the USD-free Python binding."""

import unittest

try:
    import _smlib_core as sm
except ImportError:
    sm = None


@unittest.skipIf(sm is None, "_smlib_core has not been built")
class CoreBindingTests(unittest.TestCase):
    def test_profile_is_reusable_and_returns_a_copy(self):
        vertices = [(0, 0), (0, 3), (2, 3), (2, 0), (0, 0)]
        profile = sm.create_profile(vertices)
        self.assertIsInstance(profile, sm.Profile)
        original = profile.points
        self.assertEqual(len(original), 4)
        changed = profile.points
        changed[0][0] = 100
        self.assertEqual(profile.points, original)
        first = sm.create_extrude(profile, (0, 0, 0), -1, 1)
        second = sm.create_revolve(profile, (0, 0, 0), 0, 180)
        self.assertIsNot(first, second)
        for body in (first, second):
            self.assertTrue(body.valid)
            self.assertGreater(len(sm.tessellate(body).to_numpy_arrays()["faces"]), 0)
        self.assertEqual(profile.points, original)
        self.assertEqual(vertices[0], vertices[-1])

    def test_box_and_tessellation(self):
        box = sm.create_box((0, 0, 0), 2, 3, 4)
        self.assertTrue(box.valid)
        mesh = sm.tessellate(box, chord_height_tolerance=0.05)
        self.assertIsNotNone(mesh)
        arrays = mesh.to_numpy_arrays()
        self.assertGreater(len(arrays["points"]), 0)
        self.assertGreater(len(arrays["faces"]), 0)
        self.assertEqual(len(arrays["face_normals"]), sum(
            1 for _ in _face_counts(arrays["faces"])))

    def test_boolean_consumes_second_operand(self):
        box = sm.create_box((0, 0, 0), 4, 4, 4)
        cutter = sm.create_cylinder((2, 2, -1), 0.75, 6)
        result = sm.boolean_difference(box, cutter)
        self.assertIs(result, box)
        self.assertTrue(result.valid)
        self.assertFalse(cutter.valid)
        with self.assertRaises(ValueError):
            sm.tessellate(cutter)
        self.assertGreater(len(sm.tessellate(result).to_numpy_arrays()["faces"]), 0)


def _face_counts(stream):
    offset = 0
    while offset < len(stream):
        count = stream[offset]
        yield count
        offset += count + 1


if __name__ == "__main__":
    unittest.main()
