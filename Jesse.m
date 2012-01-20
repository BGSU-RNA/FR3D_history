
Jesse,

	Run zWriteExemplarPDB to write a single pdb file of all the exemplars, arranged 20 Angstroms apart in a grid.  You can see them by running zDisplayNT('PairExemplarPDB').  They are not labeled according to what interaction, however.  To get individual files for each pair, run zWriteExemplarPDB(1).

	zDisplayExemplars