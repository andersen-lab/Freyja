import unittest

class ExampleTests(unittest.TestCase):
    def test_success(self):
        self.assertTrue(True)

    # def test_fail(self):
    #     self.assertTrue(False)

if __name__ == '__main__':
    unittest.main()