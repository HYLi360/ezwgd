import numpy as np
import pandas as pd

from ..utils.parse import blast6reader

class ConfigColl:
    def __init__(
            self,
            maxgap1: int = 40,
            maxgap2: int = 40,
            gap_penalty: int = -1,
            over_length = 0,
            over_gap: int = 3,
            min_pvalue = 0,
            coverage_ratio: float = 0.8
            ) -> None:
        self.maxgap1 = maxgap1
        self.maxgap2 = maxgap2
        self.gap_penalty = gap_penalty
        self.over_length = over_length
        self.min_pvalue = min_pvalue
        self.over_gap = over_gap
        self.coverage_ratio = coverage_ratio
        self.grading = [50, 40, 25] 

class _base_collinearity_algo:
    def __init__(self, options: ConfigColl, anchor: pd.DataFrame):
        self.mg1 = options.maxgap1
        self.mg2 = options.maxgap2
        self.gap_penalty = options.gap_penalty
        self.over_length = options.over_length
        self.min_pvalue = float(options.min_pvalue)
        self.over_gap = options.over_gap
        self.coverage_ratio = float(options.coverage_ratio)
        self.grading = options.grading

        self.anchor = anchor
        self.pvalue = 0

    def get_matrix(self):
        """Initialize the matrix for the collinearity anchor."""
        self.anchor['usedtimes1'] = 0
        self.anchor['usedtimes2'] = 0
        self.anchor['times'] = 1
        self.anchor['score1'] = self.anchor['grading']
        self.anchor['score2'] = self.anchor['grading']
        self.anchor['path1'] = self.anchor.index.to_numpy().reshape(len(self.anchor), 1).tolist()
        self.anchor['path2'] = self.anchor['path1']
        self.anchor_init = self.anchor.copy()
        self.mat_anchor = self.anchor

    def run(self):
        """Run the main collinearity processing."""
        self.get_matrix()
        self.score_matrix()
        data = []

        # Process anchor for maxPath in the positive direction
        anchor1 = self.anchor[['loc1', 'loc2', 'score1', 'path1', 'usedtimes1']].sort_values(by=['score1'], ascending=False)
        anchor1.drop(index=anchor1[anchor1['usedtimes1'] < 1].index, inplace=True)
        anchor1.columns = ['loc1', 'loc2', 'score', 'path', 'usedtimes']
        
        while (self.over_length >= self.over_gap or len(anchor1) >= self.over_gap):
            if self.max_path(anchor1):
                if self.pvalue > self.min_pvalue:
                    continue
                data.append([self.path, self.pvalue, self.score])

        # Process anchor for maxPath in the negative direction
        anchor2 = self.anchor[['loc1', 'loc2', 'score2', 'path2', 'usedtimes2']].sort_values(by=['score2'], ascending=False)
        anchor2.drop(index=anchor2[anchor2['usedtimes2'] < 1].index, inplace=True)
        anchor2.columns = ['loc1', 'loc2', 'score', 'path', 'usedtimes']

        while (self.over_length >= self.over_gap) or (len(anchor2) >= self.over_gap):
            if self.max_path(anchor2):
                if self.pvalue > self.min_pvalue:
                    continue
                data.append([self.path, self.pvalue, self.score])

        return data

    def score_matrix(self):
        """Calculate the scoring matrix for the anchor."""
        for index, row, col in self.anchor[['loc1', 'loc2']].itertuples():
            # Get anchor within a certain range
            anchor = self.anchor[(self.anchor['loc1'] > row) & 
                                 (self.anchor['loc2'] > col) & 
                                 (self.anchor['loc1'] < row + self.mg1) & 
                                 (self.anchor['loc2'] < col + self.mg2)]
            
            row_i_old, gap = row, self.mg2
            for index_ij, row_i, col_j, grading in anchor[['loc1', 'loc2', 'grading']].itertuples():
                if col_j - col > gap and row_i > row_i_old:
                    break
                score = grading + (row_i - row + col_j - col) * self.gap_penalty
                score1 = score + self.anchor.at[index, 'score1']
                if score > 0 and self.anchor.at[index_ij, 'score1'] < score1:
                    self.anchor.loc[index_ij, 'score1'] = score1
                    self.anchor.loc[index, 'usedtimes1'] += 1
                    self.anchor.loc[index_ij, 'usedtimes1'] += 1
                    self.anchor.loc[index_ij, 'path1'] = self.anchor.loc[index, 'path1'] + [index_ij]
                    gap = min(col_j - col, gap)
                    row_i_old = row_i

        # Reverse processing to handle negative direction
        anchor_reverse = self.anchor.sort_values(by=['loc1', 'loc2'], ascending=[False, True])
        for index, row, col in anchor_reverse[['loc1', 'loc2']].itertuples():
            anchor = anchor_reverse[(anchor_reverse['loc1'] < row) & 
                                    (anchor_reverse['loc2'] > col) & 
                                    (anchor_reverse['loc1'] > row - self.mg1) & 
                                    (anchor_reverse['loc2'] < col + self.mg2)]
            
            row_i_old, gap = row, self.mg2
            for index_ij, row_i, col_j, grading in anchor[['loc1', 'loc2', 'grading']].itertuples():
                if col_j - col > gap and row_i < row_i_old:
                    break
                score = grading + (row - row_i + col_j - col) * self.gap_penalty
                score2 = score + self.anchor.at[index, 'score2']
                if score > 0 and self.anchor.at[index_ij, 'score2'] < score2:
                    self.anchor.at[index_ij, 'score2'] = score2
                    self.anchor.at[index, 'usedtimes2'] += 1
                    self.anchor.at[index_ij, 'usedtimes2'] += 1
                    self.anchor.at[index_ij, 'path2'] = self.anchor.at[index, 'path2'] + [index_ij]
                    gap = min(col_j - col, gap)
                    row_i_old = row_i

    def max_path(self, anchor):
        """Find the maximum path for the given anchor."""
        if len(anchor) == 0:
            self.over_length = 0
            return False
        
        # Initialize path score and index
        self.score, self.path_index = anchor.loc[anchor.index[0], ['score', 'path']]
        self.path = anchor[anchor.index.isin(self.path_index)]
        self.over_length = len(self.path_index)
        
        # Check if the block overlaps with other blocks
        if self.over_length >= self.over_gap and len(self.path) / self.over_length > self.coverage_ratio:
            anchor.drop(index=self.path.index, inplace=True)
            [loc1_min, loc2_min], [loc1_max, loc2_max] = self.path[['loc1', 'loc2']].agg(['min', 'max']).to_numpy()

            # Calculate p-value
            gap_init = self.anchor_init[(loc1_min <= self.anchor_init['loc1']) & 
                                        (self.anchor_init['loc1'] <= loc1_max) & 
                                        (loc2_min <= self.anchor_init['loc2']) & 
                                        (self.anchor_init['loc2'] <= loc2_max)].copy()
            
            self.pvalue = self.pvalue_estimated(gap_init, loc1_max - loc1_min + 1, loc2_max - loc2_min + 1)
            self.path = self.path.sort_values(by=['loc1'], ascending=[True])[['loc1', 'loc2']]
            return True
        else:
            anchor.drop(index=anchor.index[0], inplace=True)
        return False

    def pvalue_estimated(self, gap, L1, L2):
        """Estimate p-value based on the given gap and lengths."""
        N1 = gap['times'].sum()
        N = len(gap)
        self.anchor_init.loc[gap.index, 'times'] += 1
        m = len(self.path)
        a = (1 - self.score / m / self.grading[0]) * (N1 - m + 1) / N * (L1 - m + 1) * (L2 - m + 1) / L1 / L2
        return round(a, 4)


class collinearity(_base_collinearity):
    def __init__(
            self,
            options: ConfigColl,
            anchor: pd.DataFrame,
            process: int,
            blast_reverse: bool = False
            ):
        super().__init__(options, anchor)
        self.process = process
        self.blast_reverse = blast_reverse

    def deal_blast_for_chromosomes(self, blast_path: str, rednum, repeat_number):
        # raw blast6 format: 
        # "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"
        blast = blast6reader(blast_path)
        bluenum = rednum
        def assign_grading(group):
            group['cumcount'] = group.groupby(1).cumcount()
            group = group[group['cumcount'] <= repeat_number]
            group['grading'] = pd.cut(
                group['cumcount'],
                bins=[-1, 0, bluenum, repeat_number],
                labels=self.grading,
                right=True
            )
            return group
        newblast = blast.groupby(['chr1', 'chr2']).apply(assign_grading).reset_index(drop=True)
        newblast['grading'] = newblast['grading'].astype(int)
        return newblast[newblast['grading'] > 0]
    
    def deal_blast_for_genomes(self, blast, rednum, repeat_number):
        # Initialize the grading column
        blast['grading'] = 0
        
        # Define the blue number as the sum of rednum and the predefined constant
        bluenum = 4 + rednum
        
        # Get the indices for each group by sorting the 11th column in descending order
        index = [group.sort_values(by=[11], ascending=[False])[:repeat_number].index.tolist()
                for name, group in blast.groupby([0])]
        
        # Split the indices into red, blue, and gray groups
        reddata = np.array([k[:rednum] for k in index], dtype=object)
        bluedata = np.array([k[rednum:bluenum] for k in index], dtype=object)
        graydata = np.array([k[bluenum:repeat_number] for k in index], dtype=object)
        
        # Concatenate the results into flat lists
        redindex = np.concatenate(reddata) if reddata.size else []
        blueindex = np.concatenate(bluedata) if bluedata.size else []
        grayindex = np.concatenate(graydata) if graydata.size else []

        # Update the grading column based on the group indices
        blast.loc[redindex, 'grading'] = self.grading[0]
        blast.loc[blueindex, 'grading'] = self.grading[1]
        blast.loc[grayindex, 'grading'] = self.grading[2]

        # Return only the rows with non-zero grading
        return blast[blast['grading'] > 0]

    def run(self):
        # Read and process lens files
        lens1 = base.newlens(self.lens1, 'order')
        lens2 = base.newlens(self.lens2, 'order')
        # Read and process gff files
        gff1 = base.newgff(self.gff1)
        gff2 = base.newgff(self.gff2)
        # Filter gff data based on lens indices
        gff1 = gff1[gff1['chr'].isin(lens1.index)]
        gff2 = gff2[gff2['chr'].isin(lens2.index)]
        # Process blast data

        blast = base.newblast(self.blast, int(self.score), float(self.evalue),gff1, gff2, self.blast_reverse)

        # Map positions and chromosome information
        blast['loc1'] = blast[0].map(gff1[self.position])
        blast['loc2'] = blast[1].map(gff2[self.position])
        blast['chr1'] = blast[0].map(gff1['chr'])
        blast['chr2'] = blast[1].map(gff2['chr'])
        # Apply blast filtering and grading
        if self.comparison.lower() == 'genomes':
            blast = self.deal_blast_for_genomes(blast, int(self.multiple), int(self.repeat_number))
        if self.comparison.lower() == 'chromosomes':
            blast = self.deal_blast_for_chromosomes(blast, int(self.multiple), int(self.repeat_number))
        print(f"The filtered homologous gene pairs are {len(blast)}.\n")
        if len(blast) < 1:
            print("Stopped!\n\nIt may be that the id1 and id2 in the BLAST file do not match with (gff1, lens1) and (gff2, lens2).")
            sys.exit(1)
        # Group blast data by 'chr1' and 'chr2'
        total = []
        for (chr1, chr2), group in blast.groupby(['chr1', 'chr2']):
            total.append([chr1, chr2, group])
        del blast, group
        gc.collect()
        # Determine chunk size for multiprocessing
        n = int(np.ceil(len(total) / float(self.process)))
        result, data = '', []
        try:
            # Initialize multiprocessing Pool
            pool = Pool(self.process)
            for i in range(0, len(total), n):
                # Apply single_pool function asynchronously
                data.append(pool.apply_async(
                    self.single_pool, args=(total[i:i + n], gff1, gff2, lens1, lens2)
                ))
            pool.close()
            pool.join()
        except:
            pool.terminate()
        for k in data:
            # Collect results from async tasks
            text = k.get()
            if text:
                result += text
        # Write final output to file
        result = re.split('\n', result)
        fout = open(self.savefile, 'w')
        num = 1
        for line in result:
            if re.match(r"# Alignment", line):
                # Replace alignment number
                s = f'# Alignment {num}:'
                fout.write(s + line.split(':')[1] + '\n')
                num += 1
                continue
            if len(line) > 0:
                fout.write(line + '\n')
        fout.close()
        sys.exit(0)

    def single_pool(self, group, gff1, gff2, lens1, lens2):
        text = ''
        for bk in group:
            chr1, chr2 = str(bk[0]), str(bk[1])
            print(f'Running {chr1} vs {chr2}')
            # Extract and sort points
            points = bk[2][['loc1', 'loc2', 'grading']].sort_values(
                by=['loc1', 'loc2'], ascending=[True, True]
            )
            # Initialize collinearity analysis
            collinearity = improvedcollinearity.collinearity(
                self.options, points)
            data = collinearity.run()
            if not data:
                continue
            # Extract gene information
            gf1 = gff1[gff1['chr'] == chr1].reset_index().set_index('order')[[1, 'strand']]
            gf2 = gff2[gff2['chr'] == chr2].reset_index().set_index('order')[[1, 'strand']]
            n = 1
            for block, evalue, score in data:
                if len(block) < self.over_gap:
                    continue
                # Map gene names and strands
                block['name1'] = block['loc1'].map(gf1[1])
                block['name2'] = block['loc2'].map(gf2[1])
                block['strand1'] = block['loc1'].map(gf1['strand'])
                block['strand2'] = block['loc2'].map(gf2['strand'])
                block['strand'] = np.where(
                    block['strand1'] == block['strand2'], '1', '-1'
                )
                # Prepare text output
                block['text'] = block.apply(
                    lambda x: f"{x['name1']} {x['loc1']} {x['name2']} {x['loc2']} {x['strand']}\n",
                    axis=1
                )
                # Determine alignment mark
                a, b = block['loc2'].head(2).values
                mark = 'plus' if a < b else 'minus'
                # Append alignment information
                text += f'# Alignment {n}: score={score} pvalue={evalue} N={len(block)} {chr1}&{chr2} {mark}\n'
                text += ''.join(block['text'].values)
                n += 1
        return text
