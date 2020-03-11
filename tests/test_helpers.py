from unittest import TestCase
from MCR import helpers


class TestHelpers(TestCase):
    def test_flatten_tuple(self):
        tests = [
            (("a", "b", "c"), ("a", "b", "c")),
            ((((1), "b"), 0.5), (1, "b", 0.5)),
            ((), ()),
            (True, (True,))
        ]

        for test_input, desired_output in tests:
            self.assertEqual(helpers.flatten_tuple(test_input), desired_output)


    def test_backup_dir(self):
        from pathlib import Path
        dir_to_backup = Path("../tests/backup_test_dir")
        if not dir_to_backup.exists():
            dir_to_backup.mkdir()

        backup = helpers.backup_dir(dir_to_backup)
        Path(backup).rmdir()
        self.assertIn("../tests/backup_test_dir_backup_20", backup)
        self.assertIsNone(helpers.backup_dir(dir_to_backup))


    def test_zfill_enum(self):
        test = [['a', 'b', 'c'],
                [['a', '0000'], ['b', '0001'], ['c', '0002']]]

        result = helpers.zfill_enum(test[0], padding=4)

        self.assertEqual(result, test[1])


    def test_closest_node(self):
        pass #TODO


    def test_combine_files(self):
        pass #TODO


    def test_writer(self):
        pass #TODO
