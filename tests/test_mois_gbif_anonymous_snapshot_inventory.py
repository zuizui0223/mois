import unittest
from scripts.analysis.mois_gbif_anonymous_snapshot_inventory import distribute
class TestS3Inventory(unittest.TestCase):
    def test_complete_nonoverlapping_work_shards(self):
        src=[("key"+str(i),5+i*3) for i in range(32)]
        out,sizes=distribute(src,groups=4)
        self.assertEqual(len(out),len(src))
        self.assertEqual(len({x[0] for x in out}),32)
        self.assertEqual(sum(sizes),sum(x[1] for x in src))
        self.assertTrue(all(0<=x[2]<4 for x in out))
    def test_balance(self):
        src=[("key"+str(i),100+i) for i in range(40)]
        _,sizes=distribute(src,groups=8)
        self.assertLess(max(sizes)-min(sizes),300)
    def test_bad_group_rejected(self):
        with self.assertRaises(ValueError):
            distribute([("a",100)],groups=0)
if __name__=="__main__":unittest.main()
