from ezwgd.utils.parse import blast6reader, Genome

from ezwgd.coll._corecoll import PowerColl

blast_res = blast6reader('dataset_test/utils/blast_res.f6')
genome1 = Genome('dataset_test/Pveri.fna', 'dataset_test/Pveri.gff3')
genome2 = Genome('dataset_test/Pvulg.fna', 'dataset_test/Pvulg.gff3')

PowerColl(blast_res, genome1, genome2).chunking()
