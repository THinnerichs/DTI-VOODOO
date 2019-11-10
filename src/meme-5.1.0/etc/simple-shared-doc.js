//
// simple-shared-doc.js
//
// Function to replace the innerHTML of element "id" with the HTML indicated by "doc_type".
// Easier to read and update than the more flexible approach in shared-doc.js. 
//
function print_doc(id, doc_type) {
  var html;
  switch (doc_type) {
    case 'motif-consensus':
      html = `
	<p id="consensus_doc"> 
	   A <b>consensus sequence</b> is constructed from each column in a
	   motif's frequency matrix using the <b>"50% rule"</b>
	   as follows:
	</p>
	<ol>
	  <li>The letter frequencies in the column are sorted in decreasing order.</li>
	  <li>Letters with frequency less 50% of the maximum are discarded.</li>
	  <li>The letter used in this position in the consensus sequence is determined
	  by the first rule below that applies:</li>
	  <ul>
	    <li>If there is only one letter left, or if the remaining letters exactly match
	    an ambiguous symbol in the alphabet, the <b>letter</b> or <b>ambiguous symbol</b>,
	    respectively, is used.</li>
	    <li>Otherwise, if the remaining set contains at least 50% of the core
	    symbols in the alphabet, the alphabet's <b>wildcard</b>
	    (e.g., "N" for DNA or RNA, and "X" for protein) is used.</li>
	    <li>Otherwise, the letter with the <b>maximum frequency</b> is used.</li>
	  </ul>
	</ol>
      `;
      document.getElementById(id).insertAdjacentHTML('beforeend', html);
      break;

    case 'momo-tsv-description':
      html = `
        <p>
          MoMo outputs a tab-separated values (TSV) file ('momo.tsv') that 
	  contains one line for each motif found and that is suitable 
	  for parsing by scripts and viewing with Excel.
          The first line in the file contains the (tab-separated) names of the fields.
	  Your command line is given at the end of the TSV file in a comment line 
	  starting with the character '#'.
          The names and meanings of each of the fields are described in the table below.
          Not all fields are present for all algorithms, and the field numbers are 
	  indicated for each algorithm in the first three columns of the table.
        </p>
        <table class="dark" style="width:100%" border=1>
          <tr> <th colspan=3>Field Number</th> <th rowspan=2>Field<br>Name</th> <th rowspan=2>Field<br>Contents</th> </tr>
          <tr> <th>motif-x</th> <th>MoDL</th> <th>simple</th> </tr>
          <tr> <td>1</td> <td>1</td> <td>1</td> <td>mod</td> <td>` + get_doc_text('momo-mod') + `</td> </tr>
          <tr> <td>2</td> <td>2</td> <td>2</td> <td>motif</td> <td>` + get_doc_text('momo-motif') + `</td> </tr>
          <tr> <td>3</td> <td>3</td> <td>3</td> <td>regexp</td> <td>` + get_doc_text('momo-regexp') + `</td> </tr>
          <tr> <td>4</td> <td>4</td> <td> </td> <td>score</td> <td>` + get_doc_text('momo-score') + `</td> </tr>
          <tr> <td>5</td> <td>5</td> <td>4</td> <td>fg_match</td> <td>` + get_doc_text('momo-fg') + `</td> </tr>
          <tr> <td>6</td> <td>6</td> <td> </td> <td>fg_size</td> <td>` + get_doc_text('momo-fg-size') + `</td> </tr>
          <tr> <td>7</td> <td>7</td> <td> </td> <td>bg_match</td> <td>` + get_doc_text('momo-bg') + `</td> </tr>
          <tr> <td>8</td> <td>8</td> <td> </td> <td>bg_size</td> <td>` + get_doc_text('momo-bg-size') + `</td> </tr>
          <tr> <td>9</td> <td>9</td> <td> </td> <td>fg/bg</td> <td>` + get_doc_text('momo-fold') + `</td> </tr>
          <tr> <td>10</td> <td>10</td> <td> </td> <td>unadjusted_p-value<td>` + get_doc_text('momo-unadjusted-p') + `</td> </tr>
          <tr> <td>11</td> <td> </td> <td> </td> <td>tests</td> <td>` + get_doc_text('momo-n-tests') + `</td> </tr>
          <tr> <td>12</td> <td> </td> <td> </td> <td>adjusted_p-value</td> <td>` + get_doc_text('momo-adjusted-p') + `</td> </tr>
	</table>
      `;
      document.getElementById(id).insertAdjacentHTML('beforeend', html);
      break;

    case 'momo-meme-output':
      html = `
	<p>
          MoMo outputs a file ('momo.txt') that contains each discovered motif 
          in <a href="` + site_url + `/doc/meme-format.html">MEME motif format</a>
          that is suitable for use with other MEME Suite programs.
	  MoMo creates these position-frequency matrix motifs by aligning the foreground 
	  peptides matching the motif and computing the position-frequency matrix.  
	  No pseudo-counts are added.  
	<p>
	</p>
          For motif-x motifs where the adjusted
	  <i>p</i>-value can be accurately calculated (e.g., you did not specify the
	  <code>--db-background</code> option), MoMo writes the value of 
	  <code>tests</code> times the <code>unadjusted_p-value</code> 
	  in the <i>E</i>-value field. In all other cases, MoMo writes
	  '1' in the <i>E</i>-value field of the motif.
	</p>
      `;
      document.getElementById(id).insertAdjacentHTML('beforeend', html);
      break;

    case 'ame-tsv-description':
      html = `
        <p>
          AME outputs a tab-separated values (TSV) file ('ame.tsv') that 
	  contains one line for each motif found to be significantly enriched,
	  sorted in order of decreasing statistical significance.
	  The first line in the file contains the (tab-separated) names of the fields.
	  Your command line is given at the end of the TSV file in a comment 
	  line starting with the character '#'.
          The names and meanings of each of the fields are described in the table below.
	  Not all fields are present for all types of enrichment analysis,
	  and the field numbers are indicated for each type of analysis
	  in the first six columns of the table.
	</p>
        <table class="dark" style="width:100%" border=1>
          <tr> <th colspan=6>Field Number</th> <th rowspan=2>Field<br>Name</th> <th rowspan=2>Field<br>Contents</th> </tr>
          <tr> <th>fisher</th> <th>ranksum</th> <th>3dmhg</th> <th>4dmhg</th> <th>pearson</th> <th>spearman</th> </tr> 
    	  <tr> <td colspan=6 class='ctr'>1</td> <td>rank</td> <td>The rank of the significance of the motif in the sorted results.</td> </tr>
    	  <tr> <td colspan=6 class='ctr'>2</td> <td>motif_DB</td> <td> ` + get_doc_text('motif-db', 'the motif.') + `</td> </tr>
    	  <tr> <td colspan=6 class='ctr'>2</td> <td>motif_ID</td> <td> ` + get_doc_text('motif-id') + `</td> </tr>
    	  <tr> <td colspan=6 class='ctr'>4</td> <td>motif_alt_ID</td> <td> ` + get_doc_text('motif-alt-id') + `</td> </tr>
    	  <tr> <td colspan=6 class='ctr'>5</td> <td>consensus</td> <td> ` + get_doc_text('motif-cons') + `</td> </tr>
    	  <tr> <td colspan=6 class='ctr'>6</td> <td>p-value</td> <td> ` + get_doc_text('ame-pvalue') + `</td> </tr>
    	  <tr> <td colspan=6 class='ctr'>7</td> <td>adj_p-value</td> <td> ` + get_doc_text('ame-adj-pvalue') + `</td> </tr>
    	  <tr> <td colspan=6 class='ctr'>8</td> <td>E-value</td> <td> ` + get_doc_text('ame-evalue') + `</td> </tr>
    	  <tr> <td colspan=6 class='ctr'>9</td> <td>tests</td> <td>The number of tests performed; used in correcting the <i>p</i>-value</td> </tr>
	  <tr> <td colspan=8></td> </tr>
          <tr> <td colspan=4 class='ctr'>10</td> <td></td> <td></td> <td>FASTA_max</td> <td>The optimal threshold for <b>labeling</b> sequences as positive;
            sequences with FASTA score less than or equal to the threshold are labeled as positive;
            this field will be contain the number of primary sequences if you provided
            a control file (using option <code>--control</code>); this field will contain the size of the
            partition if you specified one (using option <code>--fix-partition</code>).</td> </tr>
          <tr> <td colspan=4 class='ctr'>11</td> <td></td> <td></td> <td>pos</td> <td>The number of sequences <b>labeled</b> as positive.</td> </tr>
          <tr> <td colspan=4 class='ctr'>12</td> <td></td> <td></td> <td>neg</td> <td>The number of sequences <b>labeled</b> as negative.</td> </tr>
	  <tr> <td colspan=8></td> </tr>
          <tr> <td class='ctr'>13</td> <td></td> <td></td> <td></td> <td></td> <td></td> <td>PWM_min</td> <td>The optimal threshold on PWM score for <b>classifying</b> sequences as positive;
                sequences with PWM score greater than or equal to the threshold are classified as positive.</td> </tr>
          <tr> <td class='ctr'>14</td> <td></td> <td></td> <td></td> <td></td> <td></td> <td>TP</td> <td>The number of true positive sequences: sequences both <b>labeled</b> and <b>classified</b> as positive</td> </tr>
          <tr> <td class='ctr'>15</td> <td></td> <td></td> <td></td> <td></td> <td></td> <td>%TP</td> <td>The percentage of true positive sequences: percentage of sequences <b>labeled</b> positive and <b>classified</b> as positive.</td> </tr>
          <tr> <td class='ctr'>16</td> <td></td> <td></td> <td></td> <td></td> <td></td> <td>FP</td> <td>The number of false positive sequences: sequences <b>labeled</b> negative but <b>classified</b> as positive.</td> </tr>
          <tr> <td class='ctr'>17</td> <td></td> <td></td> <td></td> <td></td> <td></td> <td>%TP</td> <td>The percentage of false positive sequences: sequences <b>labeled</b> negative but <b>classified</b> as positive.</td> </tr>
	  <tr> <td colspan=8></td> </tr>
          <tr> <td></td> <td class='ctr'>13</td> <td></td> <td></td> <td></td> <td></td> <td>U</td> <td>The value of the Mann-Whitney <i>U</i> statistic.</td> </tr>
          <tr> <td></td> <td class='ctr'>14</td> <td></td> <td></td> <td></td> <td></td> <td>pleft</td> <td>The left-tailed <i>p</i>-value of the rank-sum test.</td> </tr>
          <tr> <td></td> <td class='ctr'>15</td> <td></td> <td></td> <td></td> <td></td> <td>pright</td> <td>The right-tailed <i>p</i>-value of the rank-sum test.</td> </tr>
          <tr> <td></td> <td class='ctr'>16</td> <td></td> <td></td> <td></td> <td></td> <td>pboth</td> <td>The two-tailed <i>p</i>-value of the rank-sum test.</td> </tr>
          <tr> <td></td> <td class='ctr'>17</td> <td></td> <td></td> <td></td> <td></td> <td>adj_pleft</td> <td>The left-tailed <i>p</i>-value of the rank-sum test, adjusted for multiple tests.</td> </tr>
          <tr> <td></td> <td class='ctr'>18</td> <td></td> <td></td> <td></td> <td></td> <td>adj_pright</td> <td>The right-tailed <i>p</i>-value of the rank-sum test, adjusted for multiple tests.</td> </tr>
          <tr> <td></td> <td class='ctr'>19</td> <td></td> <td></td> <td></td> <td></td> <td>adj_both</td> <td>The two-tailed <i>p</i>-value of the rank-sum test, adjusted for multiple tests.</td> </tr>
	  <tr> <td colspan=8></td> </tr>
	  <tr> <td></td> <td></td> <td></td> <td></td> <td class='ctr'>10</td> <td></td> <td>Pearson_CC</td> <td>The correlation coefficient of the PWM and FASTA scores of positive sequences.</td> </tr>
	  <tr> <td></td> <td></td> <td></td> <td></td> <td class='ctr'>11</td> <td></td> <td>mean_squared_error</td> <td>The mean-squared error of the regression line between PWM and FASTA scores.</td> </tr>
	  <tr> <td></td> <td></td> <td></td> <td></td> <td class='ctr'>12</td> <td></td> <td>slope</td> <td>The slope of the regression line.</td> </tr>
	  <tr> <td></td> <td></td> <td></td> <td></td> <td class='ctr'>13</td> <td></td> <td>intercept</td> <td>The y-intercept of the regression line.</td> </tr>
	  <tr> <td colspan=8></td> </tr>
	  <tr> <td></td> <td></td> <td></td> <td></td> <td></td> <td class='ctr'>10</td> <td>Spearman_CC</td> <td>The correlation coefficient of the PWM and FASTA ranks of positive sequences./td> </tr>
        </table>
	  </ul>
	</p>
      `;
      document.getElementById(id).insertAdjacentHTML('beforeend', html);
      break;

    case 'ame-sequences-tsv':
      html = `
        <p>AME outputs a tab-separated values (TSV) file ('sequences.tsv') containing one line for
	each sequence classified as 'positive' by AME for each significant motif.
	The lines are grouped by motif, and groups are separated by a line
	starting with the character "#".
        The first line in the file contains the (tab-separated) names of the fields.
        The names and meanings of each of the fields are described in the table below.
        </p>
        <table class="dark" style="width:100%" border=1>
          <tr> <th>field</th> <th>name</th> <th>contents</th> </tr>
	  <tr> <td>1</td> <td>motif_DB</td> <td>` + get_doc_text('motif-db', 'the motif.') + `</td> </tr>
	  <tr> <td>2</td> <td>motif_ID</td> <td>the ID of the motif</td> </tr>
	  <tr> <td>3</td> <td>seq_ID</td> <td>the ID of the sequence</td> </tr>
	  <tr> <td>4</td> <td><i>label_score</i> (either FASTA_score or PWM_score)</td> <td>the value of the score used to label it as positive</td> </tr>
	  <tr> <td>5</td> <td><i>class_score</i> (either PWM_score or FASTA_score)</td> <td>the value of the score used to classify it as positive</td> </tr>
	  <tr> <td>6</td> <td>class</td> <td>whether the sequence is a true positive, 'tp', or a false positive, 'fp'</td> </tr>
        </table>
      `;
      document.getElementById(id).insertAdjacentHTML('beforeend', html);
      break;

    case 'centrimo-results-tsv':
      html = `
	<p>
          CentriMo outputs a tab-separated values (TSV) file ('centrimo.tsv') that contains one line for each 
	  region found to be significantly enriched for a motif.
	  The lines are sorted in order of decreasing statistical significance.
          The first line in the file contains the (tab-separated) names of the fields.
          Your command line is given at the end of the file in a comment line starting with the
          character '#'.
          The names and meanings of each of the fields, which depend on whether or not you provide 
	  control sequences to CentriMo, are described in the table below.
	</p>
        <table class="dark" style="width:100%" border=1>
          <tr> <th>field</th> <th>name</th> <th>contents</th> </tr>
	  <tr> <td>1</td> <td>db_index</td> <td>The index of the motif file that contains the motif.  Motif
	      files are numbered in the order the appeared in the command line.</td> </tr>
	  <tr> <td>2</td> <td>motif_id</td> <td> ` + get_doc_text('motif-id') + `
	      If more than one motif has the same ID, CentriMo uses only the first such motif.
	      The name is single-quoted and preceded with '+' or '-' if you scanned separately with 
	      the reverse complement motif (using the <code>--sep</code> option).</td> </tr>
	  <tr> <td>3</td> <td>motif_alt_id</td> <td> ` + get_doc_text('motif-alt-id') + `</td> </tr>
	  <tr> <td>4</td> <td>consensus</td> <td> ` + get_doc_text('motif-cons') + `</td> </tr>
	  <tr> <td>5</td> <td>E-value</td> <td> ` + get_doc_text('centrimo-evalue') + `</td> </tr>
	  <tr> <td>6</td> <td>adj_p-value</td> <td> ` + get_doc_text('centrimo-adj-pvalue') + `
	      By default, a <i>p</i>-value is calculated by using the one-tailed binomial
	      test on the number of sequences with a match to the motif
	      that have their best match in the reported region;
              if you provided control sequences, the <i>p</i>-value of Fisher\'s exact test on the enrichment of 
              best matches in the positive sequences relative to the negative sequences is computed instead;
	      if you used the <code>--cd</code> option, the <i>p</i>-value is the probability that the average 
	      distance between the best site and the sequence center would be as low or lower than observed, 
	      computed using the cumulative Bates distribution, optimized over different score thresholds.
	      In all cases, the reported <i>p</i>-value has been adjusted for the number of regions
	      and/or score thresholds tested.</td> </tr>
	  <tr> <td>7</td> <td>log_adj_p-value</td> <td>Log of adjusted <i>p</i>-value.</td> </tr>
	  <tr> <td>8</td> <td>bin_location</td> <td>Location of the center of the most enriched region, or
		0 if you used the <code>--cd</code> option.
	  <tr> <td>9</td> <td>bin_width</td> <td> ` + get_doc_text('centrimo-bin-width') + `</td> </tr>
	  <tr> <td>10</td> <td>total_width</td> <td>The maximum number of regions possible for this motif
              <br>&nbsp;&nbsp;
	      round(sequence_length - motif_length + 1)/2,<br>
	      or the number of places the motif will fit if you used the <code>--cd</code> option.</td> </tr>
	  <tr> <td>11</td> <td>sites_in_bin</td> <td> ` + get_doc_text('centrimo-sites-in-bin') + `</td> </tr>
	  <tr> <td>12</td> <td>total_sites</td> <td>The number of sequences containing a match to the motif
	      above the score threshold.
	  <tr> <td>13</td> <td>p_success</td> <td>The probability of a random match falling into the enriched region:
	      <br>&nbsp;&nbsp;
	      bin_width / total_width</td> </tr>
	  <tr> <td>14</td> <td>p-value</td> <td>The uncorrected <i>p</i>-value before it gets adjusted for the
	      number of multiple tests to give the adjusted <i>p</i>-value.</td> </tr>
	  <tr> <td>15</td> <td>mult_tests</td> <td> ` + get_doc_text('centrimo-mult-tests') + `</td> </tr>
	  <tr> <th colspan=3>The following additional columns are present when you provide control sequences to CentriMo
	  (using the <code>--neg</code> option).</th> </tr>
	  <tr> <td>16</td> <td>neg_sites_in_bin</td> <td>The number of negative sequences where the best
	      match to the motif falls in the reported region.
	      This value is rounded but the underlying value may contain fractional counts.
	      Note: This number may be less than the number of negative have a best match in the region.
	      The reason for this is that a sequence may have many matches that score equally best.
	      If n matches have the best score in a sequence, 1/n is added to the
	      appropriate bin for each match.</td> </tr>
	  <tr> <td>17</td> <td>neg_sites</td> <td>The number of negative sequences containing a match to the
	      motif above the minimum score threshold.
	      When score optimization is enabled the score threshold may be raised
	      higher than the minimum.</td> </tr>
	  <tr> <td>18</td> <td>neg_adj_pvalue</td> <td>The probability that any tested region in the negative
	      sequences would be as enriched for best matches to this motif
	      according to the Binomial test.</td> </tr>
	  <tr> <td>19</td> <td>log_neg_adj_pvalue</td> <td>Log of negative adjusted <i>p</i>-value.</td> </tr>
	  <tr> <td>20</td> <td>fisher_adj_pvalue</td> <td>Fisher adjusted <i>p</i>-value before it gets adjusted for the
	      number of motifs in the input files(s).</td> </tr>
	  <tr> <td>21</td> <td>log_fisher_adj_pvalue</td> <td>Log of Fisher adjusted <i>p</i>-value.</td> </tr>
        </table>
      `;
      document.getElementById(id).insertAdjacentHTML('beforeend', html);
      break;

    case 'centrimo-sites-txt':
      html = `
	<p>
          CentriMo outputs a text file ('site_counts.txt') that contains,
	  for each motif, pairs of values (bin_position, site_count),
	  or triples of values (bin_position, site_count, neg_site_count) if
	  you provided control sequences to CentriMo.
	  This data can be used to plot the density of motif best matches (sites)
	  along the input sequences.  Fractional counts are possible if multiple (n)
	  bins contain the best match for a given sequence, with each bin  
	  receiving an incremental count of 1/n. 
	</p>
	<p>
	  The data for each motif begins with a header line with the format:
	  <br>&nbsp&nbsp
		DB &lt;db_number&gt; MOTIF &lt;id&gt; &lt;alt&gt;
	  </br>
	  where &lt;id&gt; and &lt;alt&gt; are as described above.
	  The following lines (up to the next header line) 
	  each contain a single value-pair or value-triple for the motif
	  named in the header line.
	</p>
      `;
      document.getElementById(id).insertAdjacentHTML('beforeend', html);
      break;

    case 'meme-chip-results-tsv':
      html = `
	<p>
          MEME-ChIP outputs a tab-separated values (TSV) file ('summary.tsv') that 
	  contains one line for each motif found by MEME-ChIP.
	  The lines are sorted in order of decreasing statistical significance.
	  The first line in the file contains the (tab-separated) names of the fields.
	  Your command line is given at the end of the file in a comment line starting with the
	  character '#'.
          The names and meanings of the fields in the are described in the table below. 
	</p>
        <table class="dark" style="width:100%" border=1>
          <tr> <th>field</th> <th>name</th> <th>contents</th> </tr>
          <tr> <td>1</td> <td>MOTIF_INDEX</td> <td>The index of the motif in the "Motifs in MEME text format" file ('combined.meme') 
		output by MEME-ChIP.</td> </tr>
          <tr> <td>2</td> <td>MOTIF_SOURCE</td> <td>The name of the program that found the <i>de novo</i> motif, or the
                name of the motif file containing the known motif.</td> </tr>
          <tr> <td>3</td> <td>MOTIF_ID</td> <td> ` + get_doc_text('motif-id') + `</td> </tr>
          <tr> <td>4</td> <td>ALT_ID</td> <td> ` + get_doc_text('motif-alt-id') + `</td> </tr>
          <tr> <td>5</td> <td>CONSENSUS</td> <td>The ID of the <i>de novo</i> motif, or a consensus sequence
                computed from the letter frequencies in the known motif
                (as described <a href="#consensus_doc">below</a>).</td> </tr>
          <tr> <td>6</td> <td>WIDTH</td> <td>The width of the motif.</td> </tr>
          <tr> <td>7</td> <td>SITES</td> <td>The number of sites reported by the <i>de novo</i> program, or the number
                of "Total Matches" reported by CentriMo.</td> </tr>
          <tr> <td>8</td> <td>E-VALUE</td> <td>The statistical significance of the motif.</td> </tr>
          <tr> <td>9</td> <td>E-VALUE_SOURCE</td> <td>The program that reported the <i>E</i>-value.</td> </tr>
          <tr> <td>10</td> <td>MOST_SIMILAR_MOTIF</td> <td>The known motif most similar to this motif according to Tomtom.</td> </tr>
          <tr> <td>11</td> <td>URL</td> <td>A link to a description of the most similar motif, or to the known motif.</td> </tr>
        </table>
      `;
      document.getElementById(id).insertAdjacentHTML('beforeend', html);
      break;

    case 'meme-chip-combined-motifs':
      html = `
        <p>
          MEME-ChIP outputs a text file ('combined.meme') containing all the significant motifs found by MEME-ChIP.
          The motifs are in <a href="` + site_url + `/doc/meme-format.html">Minimal MEME Motif format</a>, 
	  and their IDs correspond to the motif indices given in the "Summary in TSV Format" file ('summary.tsv').
        </p>
        </p>
          <b>Note:</b> The 'nsites=' and 'E=' fields in the motif headers are only
          relevant for the MEME and DREME motifs.  For known motifs, those values do
          not refer to the number of sites in the input sequences.
        </p>
      `;
      document.getElementById(id).insertAdjacentHTML('beforeend', html);
      break;

    case 'spamo-results-tsv':
      html = `
	<p>
          SpaMo outputs a tab-separated values (TSV) file ('spamo.tsv') that contains one line for each motif
	  found to be significantly enriched.
	  The lines are grouped by secondary motif and sorted in order of decreasing statistical significance.
          The first line in the file contains the (tab-separated) names of the fields.
	  Your command line is given at the end of the file in a comment line starting with the character '#'.
	  The names and meanings of each of the fields are described in the table below. 
	</p>
	<table class="dark" style="width:100%" border=1>
	  <tr> <th>field</th> <th>name</th> <th>contents</th> </tr>
          <tr> <td>1</td> <td>prim_db</td> <td> ` + get_doc_text('motif-db', 'the primary motif.') + `</td> </tr>
          <tr> <td>2</td> <td>prim_id</td> <td> ` + get_doc_text('motif-id', 'primary') + `</td> </tr>
          <tr> <td>3</td> <td>prim_alt</td> <td> ` + get_doc_text('motif-alt-id', 'primary') + `</td> </tr>
          <tr> <td>4</td> <td>prim_cons</td> <td> ` + get_doc_text('motif-cons', 'primary') + `</td> </tr>
          <tr> <td>5</td> <td>sec_db</td> <td> ` + get_doc_text('motif-db', 'the secondary motif.') + `</td> </tr>
          <tr> <td>6</td> <td>sec_id</td> <td> ` + get_doc_text('motif-id', 'secondary') + `</td> </tr>
          <tr> <td>7</td> <td>sec_alt</td> <td> ` + get_doc_text('motif-alt-id', 'secondary') + `</td> </tr>
          <tr> <td>8</td> <td>sec_cons</td> <td> ` + get_doc_text('motif-cons', 'secondary') + `</td> </tr>
          <tr> <td>9</td> <td>trim_left</td> <td>Number of positions trimmed from left of secondary motif.</td> </tr>
          <tr> <td>10</td> <td>trim_right</td> <td>Number of positions trimmed from right of secondary motif.</td> </tr> 
          <tr> <th colspan=3>If the next three fields are not blank, the motif is redundant with a more significant ('parent') motif.</th> </tr>
          <tr> <td>11</td> <td>red_db</td> <td> ` + get_doc_text('motif-db', 'the parent motif.') + `</td> </tr>
          <tr> <td>12</td> <td>red_id</td> <td> ` + get_doc_text('motif-id', 'parent') + `</td> </tr>
          <tr> <td>13</td> <td>red_alt</td> <td> ` + get_doc_text('motif-alt-id', 'parent') + `</td> </tr>
          <tr> <td>14</td> <td>E-value</td> <td>The expected number motifs that would have least one spacing as enriched as the best spacing for this secondary. 
	    The <i>E</i>-value is the best spacing <i>p</i>-value multiplied by the number of motifs in the input database(s).</td> </tr>
          <tr> <td>15</td> <td>gap</td> <td>The distance between the edge of the primary and the (trimmed) secondary motif.</td> </tr>
          <tr> <td>16</td> <td>orient</td> <td>The (combination) of quadrants for which occurrences of this spacing are combined.</td> </tr>
          <tr> <td>17</td> <td>count</td> <td>The number of occurrences of the secondary motif with the given spacing and orientation to the primary motif.</td> </tr>
          <tr> <td>18</td> <td>total</td> <td>The total number of occurrences of the secondary motif within the margins around the best primary motif occurrence.</td> </tr>
          <tr> <td>19</td> <td>adj_p-value</td> <td>The <i>p</i>-value of the gap and orientation, adjusted for nine combinations of quadrants times the number of gaps tested (as controlled by the <code>-range</code> option).</td> </tr>
          <tr> <td>20</td> <td>p-value</td> <td>The <i>p</i>-value of the gap and orientation adjusted only for the number of gaps tested.</td> </tr>
        </table>
        <br>
      `;
      document.getElementById(id).insertAdjacentHTML('beforeend', html);
      break;

    case 'spamo-dumpseqs-tsv':
      html = `
	<p>
	  By specifying the options <code>--dumpseqs</code> or <code>--dumpsigs</code> 
	  you can have SpaMo create TSV (tab-separated values) files
	  describing the motif matches SpaMo used to make the histograms in its HTML output.
	  The files are named 
          '<code>seqs_&lt;primary_motif&gt;_&lt;secondary_db&gt;_&lt;secondary_motif&gt;.txt</code>'.
	  The rows in each file are sorted by sequence name.  
	  The first line in the file contains the (tab-separated) names of the fields.
	  The names and meanings of each of the fields are described in the table below. 
        </p>
	<table class="dark" style="width:100%">
	  <tr><th>field</th><th>name</th><th>contents</th></tr>
	  <tr><td>1</td><td>matches</td><td>Trimmed lowercase sequence centered on primary match with matches in uppercase.</td></tr>
	  <tr><td>2</td><td>sec_pos</td><td>Position of the secondary match within the whole sequence.</td></tr>
	  <tr><td>3</td><td>pri_match</td><td>Sequence fragment that the primary matched.</td></tr>
	  <tr><td>4</td><td>pri_strand</td><td>Strand of the primary match (+/-).</td></tr>
	  <tr><td>5</td><td>sec_match</td><td>Sequence fragment that the secondary matched.</td></tr>
	  <tr><td>6</td><td>sec_strand</td><td>Strand of the secondary match (+/-).</td></tr>
	  <tr><td>7</td><td>same_opp</td><td>The primary match on the same (s) or opposite (o) strand as the secondary.</td></tr>
	  <tr><td>8</td><td>down_up</td><td>The secondary match is downstream (d) or upstream (u) of the primary.</td></tr>
	  <tr><td>9</td><td>gap</td><td>The gap between the primary and secondary matches.</td></tr>
	  <tr><td>10</td><td>seq_name</td><td>The name of the sequence.</td></tr>
	  <tr><td>11</td><td>adj_p-value</td><td>The <i>p</i>-value of the bin containing the match, adjusted for the number of bins.</td></tr>
	  <tr><th colspan="3">If the sequence names are in UCSC Genome Browser position
	  format (e.g., "chr5:36715616-36715623"), the following additional fields will be present:</th></tr>
	  <tr><td>12</td><td>pri_bed_chr</td><td>Position of primary match in BED coordinates.</td></tr>
	  <tr><td>13</td><td>pri_bed_start</td><td>"</td></tr>
	  <tr><td>14</td><td>pri_bed_end</td><td>"</td></tr>
	  <tr><td>15</td><td>pri_browser</td><td>Position of primary match in UCSC Genome Browser coordinates.</td></tr>
	  <tr><td>16</td><td>sec_bed_chr</td><td>Position of secondary match in BED coordinates.</td></tr>
	  <tr><td>16</td><td>sec_bed_start</td><td>"</td></tr>
	  <tr><td>16</td><td>sec_bed_end</td><td>"</td></tr>
	  <tr><td>19</td><td>sec_browser</td><td>Position of secondary match in UCSC Genome Browser coordinates.</td></tr>
	</table>
      `;
      document.getElementById(id).insertAdjacentHTML('beforeend', html);
      break;

    case 'fimo-results-tsv':
      html = ` 
        <p>
          FIMO outputs a tab-separated values (TSV) file ('fimo.tsv') that contains one line for each
          significant match to a motif.
          The lines are sorted in order of decreasing statistical significance (increasing <i>p</i>-value).
          The first line contains the (tab-separated) names of the fields.
          Your command line is given at the end of the file in a comment line starting with the
          character '#'.  The names and meanings of each of the fields are described in the table below. 
	</p>
	<table class="dark" style="width:100%" border=1>
	  <tr> <th>field</th> <th>name</th> <th>contents</th> </tr>
          <tr> <td>1</td> <td>motif_id</td> <td> ` + get_doc_text('motif-id') + `</td> </tr>
          <tr> <td>2</td> <td>motif_alt_id</td> <td> ` + get_doc_text('motif-alt-id') + `</td> </tr>
          <tr> <td>3</td> <td>sequence_name</td> <td> ` + get_doc_text('sequence-id') + ` --OR-- the `
	    + get_doc_text('sequence-name') + `</td> </tr>
          <tr> <td>4</td> <td>start</td> <td> ` + get_doc_text('match-start-seq', 'motif occurrence') + ` --OR-- `
            + get_doc_text('match-start-genomic', 'motif occurrence') 
            + get_doc_text('parse-genomic-coord', 'The latter case occurs when FIMO') + `</td> </tr>
          <tr> <td>5</td> <td>stop</td> <td> ` + get_doc_text('match-stop-seq', 'motif occurrence') + ` --OR-- `
            + get_doc_text('match-stop-genomic', 'motif occurrence') 
            + get_doc_text('parse-genomic-coord', 'The latter case occurs when FIMO') + `</td> </tr>
          <tr> <td>6</td> <td>strand</td> <td>The strand '<code>+</code>' indicates the motif matched the forward
	    strand, '<code>-</code>' the reverse strand, and '<code>.</code>'
	    indicates strand is not applicable (as for amino acid sequences).</td> </tr>
          <tr> <td>7</td> <td>score</td> <td>The score for the motif occurrence. `
            + get_doc_text('motif-match-score') + `</td> </tr>
          <tr> <td>8</td> <td>p-value</td> <td> The <i>p</i>-value of the motif occurrence. `
            + get_doc_text('motif-match-p-value') + `</td> </tr>
          <tr> <td>9</td> <td>q-value</td> <td>The q-value of the motif occurrence. `
            + get_doc_text('bh-q-value', 'FIMO') + ` <b>Note:</b> This column is empty if 
	   you used the <code>--text</code> switch.</td> </tr>
          <tr> <td>10</td> <td>sequence</td> <td>The region of the sequence matched to the motif.</td> </tr>
        </table>
      `;
      document.getElementById(id).insertAdjacentHTML('beforeend', html);
      break;

    case 'fimo-results-gff3':
      html = `
        <p>
          FIMO outputs a GFF3</a> file ('fimo.gff') that contains one line for each
          significant match to a motif.
        </p>
        <p>
	  The GFF3 format is described <a href="http://gmod.org/wiki/GFF3">here</a>.
	  The 'score' reported in the FIMO GFF3 output</a> (in column 6) is<br/>
	  &nbsp;&nbsp;&nbsp;&nbsp;<code>min(1000, -10*(log10(pvalue)))</code>, <br/>
	  and the 'display name' ('Name=' tag in column 9) is composed of the contents of three 
	  fields as follows <br/>
	  &nbsp;&nbsp;&nbsp;&nbsp;<code>&lt;motif_id&gt;_&lt;sequence_name&gt;&lt;strand&gt;</code>.
        </p>
      `;
      document.getElementById(id).insertAdjacentHTML('beforeend', html);
      break;

    case 'gomo-results-tsv':
      html = `
        <p>
          GOMo outputs a tab-separated values (TSV) file ('gomo.tsv') that contains one line for each 
	  motif-GO-term pair found to be significantly enriched.
          The lines are grouped by motif and sorted in order of decreasing statistical significance.
          The first line contains the (tab-separated) names of the fields.
          Your command line is given at the end of the file in a comment line starting with the character '#'.
          The names and meanings of each of the fields are described in the table below.
        </p>
        <table class="dark" style="width:100%" border=1>
          <tr> <th>field</th> <th>name</th> <th>contents</th> </tr>
          <tr> <td>1</td> <td>Motif_Identifier</td> <td> ` + get_doc_text('motif-id') + ` </td> </tr> 
          <tr> <td>2</td> <td>GO_Term_Identifier</td> <td> ` + get_doc_text('gomo-go-term') + ` </td> </tr>
          <tr> <td>3</td> <td>GOMo_Score</td> <td> ` + get_doc_text('gomo-score') + ` </td> </tr>
          <tr> <td>4</td> <td>p-value</td> <td> ` + get_doc_text('gomo-p-value') + ` </td> </tr>
          <tr> <td>5</td> <td>q-value</td> <td> ` + get_doc_text('bh-q-value', 'GOMo') + ` </td> </tr>
        </table>
      `;
      document.getElementById(id).insertAdjacentHTML('beforeend', html);
      break;

    case 'gomo-results-xml':
      html = `
	<p>
	GOMo outputs an XML file ('gomo.xml') with the following format.
	</p>
	<table class="bordertable" border="1">
	  <tr>
	    <th>Tag</th><th>Child of</th><th>Description</th>
	  </tr>
	  <tr>
	    <td >&lt;gomo&gt;</td><td >Nothing</td>
	    <td>
	      Information about this run of GOMo.
	      <ul>
		<li>version - The version of GOMo that generated the XML file.</li>
		<li>release - The release date of the version that generated the XML.</li>
	      </ul>
	    </td>
	  </tr>
	  <tr>
	    <td >&lt;program&gt;</td>
	    <td >&lt;gomo&gt;</td>
	    <td>
	      Information about the state of the program when it ran.<br />
	      <ul>
		<li>name - name of the program.</li>
		<li>cmd - the command line passed to the program.</li>
		<li>gene_url - the url used to lookup further information on the gene ids. 
		The url has ampersands (&amp;) converted into <b>&amp;amp;</b> and the place where
		  the gene ID should be replaced by <b>!!GENEID!!</b> .</li>
		<li>outdir - the output directory that the program wrote to.</li>
		<li>clobber - true if GOMo was allowed to overwrite the output directory.</li>
		<li>text_only - true if GOMo wrote to stdout, in which case this file would
		  not exist so it must be false.</li>
		<li>use_e_values - true if GOMo used <i>E</i>-values (converted from <i>p</i>-values) as 
		  input scores, false if GOMo used gene scores.</li>
		<li>score_e_thresh - if GOMo used <i>E</i>-values then this is the threshold that 
		  GOMo assumed the worst <i>E</i>-value (<i>p</i>-value = 1.0) for the gene to smooth out noise.</li>
		<li>min_gene_count - the minimum number of genes that a GO term was annotated 
		  with before GOMo would calculate a score for it.</li>
		<li>motifs - if present then a space delimited list of the motifs that GOMo
		  calculated a score for, otherwise GOMo scored all motifs.</li>
		<li>shuffle_scores - the number of times GOMo generated a shuffled mapping of
		  gene id to gene id to be used to generate scores from the null model.</li>
		<li>q_threshold - GOMo filtered the results to only show those with a better
		(smaller) q-value.</li>
	      </ul>
	    </td>
	  </tr>
	  <tr>
	    <td>&lt;gomapfile&gt;</td>
	    <td>&lt;program&gt;</td>
	    <td>
	      Information about the GO mapping file.<br />
	      <ul>
		<li>path - the path to the mapping file.</li>
	      </ul>
	    </td>
	  </tr>
	  <tr>
	    <td>&lt;seqscorefile&gt;</td>
	    <td>&lt;program&gt;</td>
	    <td>
	      Information about the sequence scoring file.<br />
	      <ul>
		<li>path - the path to the sequence scoring file.</li>
	      </ul>
	    </td>
	  </tr>
	  <tr>
	    <td>&lt;motif&gt;</td>
	    <td>&lt;gomo&gt;</td>
	    <td>
	      Information about the motif.<br />
	      <ul>
		<li>id - the motif identifier.</li>
		<li>genecount - the number of scored sequences that were used to compute the result.</li>
	      </ul>
	    </td>
	  </tr>
	  <tr>
	    <td>&lt;goterm&gt;</td>
	    <td>&lt;motif&gt;</td>
	    <td>
	      Information about the GO term.<br />
	      <ul>
		<li>id - the GO identifier.</li>
		<li>score - the geometric mean across all species of the rank-sum test <i>p</i>-value.</li>
		<li>pvalue - the empirically calculated <i>p</i>-value.</li>
		<li>qvalue - the empirically calculated q-value.</li>
		<li>annotated - the number of genes annotated with the go term.</li>
		<li>group - the subgroup that the term belongs to. For the Gene Ontology 
		    b = biological process, c = cellular component and m = molecular function.</li>
		<li>nabove - the number of more general terms that link to this one.</li>
		<li>nbelow - the number of more specific terms that link from this one.</li>
		<li>implied - is the go term implied by other significant go terms? 
		  Allows values 'y', 'n' or 'u' (default) for yes, no or unknown.</li>
		<li>description - the GO term description.</li>
	      </ul>
	    </td>
	  </tr>
	  <tr>
	    <td>&lt;gene&gt;</td>
	    <td>&lt;goterm&gt;</td>
	    <td>
	      Information about the GO term's annotated genes for the primary species.<br />
	      <ul>
		<li>id - the gene identifier.</li>
		<li>rank - the rank of the scored gene.</li>
	      </ul>
	    </td>
	  </tr>
        </table>
      `;
      document.getElementById(id).insertAdjacentHTML('beforeend', html);
      break;

    case 'mcast-results-tsv':
      html = ` 
        <p>
          MCAST outputs a tab-separated values (TSV) file ('mcast.tsv') that contains one line for each
          significant match of a cluster of motifs to a sequence region.
          The lines are sorted in order of decreasing statistical significance (increasing <i>p</i>-value).
          The first line contains the (tab-separated) names of the fields.
          Your command line is given at the end of the file in a comment line starting with the
          character '#'.  The names and meanings of each of the fields are described in the table below.
	</p>
	<table class="dark" style="width:100%" border=1>
	  <tr> <th>field</th> <th>name</th> <th>contents</th> </tr>
          <tr> <td>1</td> <td>pattern_name</td> <td>A unique name that MCAST generates for the match. </td> </tr>
          <tr> <td>2</td> <td>sequence_name</td> <td> ` + get_doc_text('sequence-id') + ` --OR-- the `
	    + get_doc_text('sequence-name') + `</td> </tr>
          <tr> <td>3</td> <td>start</td> <td> ` + get_doc_text('match-start-seq', 'matched sequence region') + ` --OR-- `
            + get_doc_text('match-start-genomic', 'motif occurrence')
            + get_doc_text('parse-genomic-coord', 'The latter case occurs when MCAST') + `</td> </tr>
          <tr> <td>4</td> <td>stop</td> <td> ` + get_doc_text('match-stop-seq', 'matched sequence region') + ` --OR-- `
            + get_doc_text('match-stop-genomic', 'motif occurrence')
            + get_doc_text('parse-genomic-coord', 'The latter case occurs when MCAST') + `</td> </tr>
          <tr> <td>5</td> <td>score</td> <td> ` + get_doc_text('mcast-cluster-score') + ` </td> </tr>
          <tr> <td>6</td> <td>p-value</td> <td> ` + get_doc_text('mcast-cluster-p-value') + ` </td> </tr>
          <tr> <td>7</td> <td>E-value</td> <td> ` + get_doc_text('mcast-cluster-E-value') + ` </td> </tr>
          <tr> <td>8</td> <td>q-value</td> <td> ` + get_doc_text('bh-q-value', 'MCAST') + ` 
	    <b>Note:</b> This column is empty if you used the <code>--text</code> switch.</td> </tr>
          <tr> <td>9</td> <td>sequence</td> <td>The region of the sequence matched to the motif cluster.</td> </tr>
        </table>
      `;
      document.getElementById(id).insertAdjacentHTML('beforeend', html);
      break;

    case 'mcast-results-gff3':
      html = `
        <p>
          MCAST outputs a GFF3</a> file ('mcast.gff') that contains one line for each
          significant match to a motif cluster.
        </p>
        <p>
	  The GFF3 format is described <a href="http://gmod.org/wiki/GFF3">here</a>.
	  The 'score' reported in the MCAST GFF3 output</a> (in column 6) is<br/>
	  &nbsp;&nbsp;&nbsp;&nbsp;<code>min(1000, -10*(log10(pvalue)))</code>, <br/>
	  and the 'unique identifier' ('ID=' tag in column 9) is the value of the
	  &lt;pattern_name&gt; field in the MCAST results TSV format.  Following the
	  unique identifier in column 9, the <i>p</i>-value, <i>E</i>-value and q-value
	  of the match are given.
        </p>
      `;
      document.getElementById(id).insertAdjacentHTML('beforeend', html);
      break;

    case 'tomtom-results-tsv':
      html = `
	<p>
          Tomtom outputs a tab-separated values (TSV) file ('tomtom.tsv') that contains one line for each motif
	  found to be significantly enriched.
	  The lines are grouped by query motif and sorted in order of decreasing statistical significance.
          The first line contains the (tab-separated) names of the fields.
	  Your command line is given at the end of the file in a comment line starting with the character '#'.
          The names and meanings of each of the fields are described in the table below.
	</p>
	<table class="dark" style="width:100%" border=1>
	  <tr> <th>field</th> <th>name</th> <th>contents</th> </tr>
          <tr> <td>1</td> <td>Query_ID</td> <td> ` + get_doc_text('motif-id', 'query') + `</td> </tr>
          <tr> <td>2</td> <td>Target_ID</td> <td> ` + get_doc_text('motif-id', 'target') + `</td> </tr>
          <tr> <td>3</td> <td>Optimal_offset</td> <td> ` + get_doc_text('tomtom-offset') + `</td> </tr>
          <tr> <td>4</td> <td>p-value</td> <td> ` + get_doc_text('tomtom-p-value') + `</td> </tr>
          <tr> <td>5</td> <td>E-value</td> <td> ` + get_doc_text('tomtom-E-value') + `</td> </tr>
          <tr> <td>6</td> <td>q-value</td> <td> ` + get_doc_text('bh-q-value', 'Tomtom') + `</td> </tr>
          <tr> <td>7</td> <td>Overlap</td> <td> ` + get_doc_text('tomtom-overlap') + `</td> </tr>
          <tr> <td>8</td> <td>Query_consensus</td> <td>A consensus sequence
                computed from the letter frequencies in the query motif
                (as described <a href="#consensus_doc">below</a>).</td> </tr>
          <tr> <td>9</td> <td>Target_consensus</td> <td>A consensus sequence
                computed from the letter frequencies in the target motif
                (as described <a href="#consensus_doc">below</a>).</td> </tr>
          <tr> <td>10</td> <td>Orientation</td> <td> ` + get_doc_text('tomtom-orientation', "<br>A value of '+' means that the target motif is as it appears in the database. A value of '-' means that the reverse complement of the target motif is shown.") + `</td> </tr>
        </table>
      `;
      document.getElementById(id).insertAdjacentHTML('beforeend', html);
      break;

    default:
      html = "Error--Unrecognized doc_type: " + doc_type;
      document.getElementById(id).insertAdjacentHTML('beforeend', html);
  }
} // print_doc

//
// Function to return the HTML text of a given type.
// This function can be used directly to document the output format (xx-output-format.html)
// and indirectly via print_doc_para for help pop-ups in the actual output HTML,
// to prevent duplication of documentation.
//
function get_doc_text(doc_type, extra) {
  var html;
  if (extra == undefined) {extra = ""};

  switch (doc_type) {

    // shared
    case 'motif-db':
      return(`
	The name of a file of motifs ("motif database file") that contains ` + extra + `
      `);
    case 'motif-id':
      return(`
	The name of the ` + extra + ` motif, which is unique in the motif database file.
      `);
    case 'motif-alt-id':
      return(`
	An alternate name for the ` + extra + ` motif, which may be provided in the motif database file.
      `);
    case 'motif-width':
      return(`
	The width of the motif. No gaps are allowed in motifs supplied to ` + extra + `
        as it only works for motifs of a fixed width.
      `);
    case 'motif-cons':
      return(`
	A consensus sequence computed from the ` + extra + ` motif (as described <a href="#consensus_doc">below</a>).
      `);
    case 'motif-match-score':
     return(`
	The score for the match of a position in a sequence to a motif is
	computed by summing the appropriate entry from each column of the
	position-dependent scoring matrix that represents the motif. ` + extra + `
     `);
    case 'motif-match-p-value':
      return(`
	The <i>p</i>-value of a motif match is the probability of a single random
	subsequence of the length of the motif
	<a href="javascript:help_refine('pop_motif_match_score')">scoring</a>
	at least as well as the observed match.
      `);
    case 'bh-q-value-method':
      return(`
        ` + extra + ` estimates q-values from all the match <i>p</i>-values 
	using the method proposed by Benjamini & Hochberg (<i>Journal of the Royal Statistical Society B</i>, 57:289-300, 1995).
	See also Storey JD, Tibshirani R. Statistical significance for
	genome-wide studies, <i>Proc. Natl. Acad. Sci. USA</i> (2003) <b>100</b>:9440&ndash;9445.
      `);
    case 'bh-q-value':
      return(`
	The minimum False Discovery Rate (FDR) required to consider this match significant.</br>` + get_doc_text('bh-q-value-method', extra) + `
      `);
    case 'sdb-name':
      return(`
	The name of the (FASTA) sequence database file.
      `);
    case 'sdb-psp':
      return(`
	The name of the position specific priors (PSP) file.
      `);
    case 'sdb-dist':
      return(`
	The name of the binned distribution of priors file.
      `);
    case 'sdb-count':
      return(`
	The number of sequences in the database.
      `);
    case 'sdb-letters':
      return(`
	The number of letters in the sequence database.
      `);
    case 'lastmod':
      return(`
	The date of the last modification to the ` + extra + ` database.
      `);
    case 'sequence-id':
      return(`
        The identifier of the sequence (from the FASTA sequence header line)` + extra + `
      `);
    case 'sequence-name':
      return(`
	` + extra + `name of the sequence extracted from the sequence identifier (in the FASTA sequence header line).<br>
	When you use the <code>--parse-genomic--coord</code> option, the sequence name ends at the
	first colon ':' (if any) present in the sequence\'s FASTA identifier.  Typically this is the
	chromosome or contig name.  With the <code>--parse-genomic--coord</code> option,
	the start and stop positions are in 0-based coordinates relative to the sequence start given 
	in the FASTA sequence identifier (just after the sequence name).</td> </tr>
      `);
    case 'sequence-desc':
      return(`
        The description appearing after the identifier of the sequence in the FASTA header line.
      `);
    case 'sequence-name':
    case 'alph-name':
      return(`
	The name of the alphabet symbol.
      `);
    case 'alph-bg':
      return(`
	The frequency of the alphabet symbol as defined by the background model.
      `);
    case 'match-start-seq':
      return(`
	The start position of the ` + extra + `; 1-based sequence coordinates.
      `);
    case 'match-stop-seq':
      return(`
	The end position of the ` + extra + `; 1-based sequence coordinates.
      `);
    case 'match-start-genomic':
      return(`
	The start position of the ` + extra + `; genomic coordinates.
      `);
    case 'match-stop-genomic':
      return(`
	The end position of the ` + extra + `; genomic coordinates.
      `);
    case 'parse-genomic-coord':
      return(`
	` + extra + ` was run with the <code>--parse-genomic-coord</code> option
	and has split the sequence identifier into sequence name, sequence start and sequence end 
	in genomic coordinates.
      `);

    // MoMo output fields
    case 'momo-logo':
      return(`
        MoMo creates a sequence logo for each motif it discovers using
        the <a href="` + site_url + `/doc/ceqlogo.html">ceqlogo</a> utility.
      `);
    case 'momo-mod':
      return(`
	The post-translationally modified residue located
        at the center of the motif.
      `);
    case 'momo-motif':
      return(`
        (A string similar to) a regular expression describing the motif.
        Lower case 'x' represents a match to any residue the non-central 
	positions.
        The central peptide and its modification weight (if any) is surrounded
        by '_' characters. If option <code>--single-motif-per-mass</code> 
	was specified, the central peptide will be represented by an uppercase 'X', 
	which also represents a match to any residue.  With the simple algorithm, 
	all non-central positions will contain an 'x', regardless of which residues appeared in
        the input sequences.  (See the Regular Expression column of the output
        for further information.) Motifs are reported in the order 
	in which they are found by MoMo, and are not sorted in any way.
      `);
    case 'momo-regexp':
      return(`
	A PERL regular expression suitable for searching for the motif in
        protein sequences using 
	<a href="` + site_url + `/doc/fasta-grep.html">fasta-grep</a> 
	or PERL scripts. The PERL wildcard character '.' matches any
	resudue. Groups of residues in square brackets '[ ]' match
	any of those residues. For the motif-x and MoDL algorithms,
	this is the regular expression for the motif as output by that
	algorithm.  For the simple algorithm,
	this regular expression includes all of the residues 
	in the motif occurrences that match the motif.
      `);
    case 'momo-score':
      return(` 
	The algorithm-dependent score of the motif.
        For algorithm motif-x, this is the sum over the significant 
	position/residue pairs of -log(p<sub>binomial</sub>);
        for algorithm MoDL, this is the increase in description length in 
	bits if the motif were to be removed from the final set of motifs.
      `);
    case 'momo-fg':
      return(`
	fg_match is the number of foreground peptides that match the motif.
        <b>Note: </b>For motif-x, the foreground counts are based on
        <b>all</b> the peptides unless the 
	<code>--harvard</code> option is specified, in which 
	case the counts are based on the <b>remaining</b> peptides.
      `);
    case 'momo-fg-size':
      return(`
	fg_size is the total number of foreground peptides with the given 
	central modification.
      `);
    case 'momo-bg':
      return(`
	bg_match is the number of background peptides that match the motif.
        <b>Note: </b>For motif-x, the background counts are based on
        <b>all</b> the peptides unless the 
	<code>--harvard</code> option is specified, in which 
	case the counts are based on the <b>remaining</b> peptides.
      `);
    case 'momo-bg-size':
      return(`
	bg_size is the total number of background peptides with the given 
	central modification.
      `);
    case 'momo-fold':
      return(`
	The fold enrichment of the foreground matches vs. the background
	matches.  This is equal to: <br>
	&nbsp;&nbsp;(fg_match / fg_size) / (bg_match / bg_size).
      `);
    case 'momo-unadjusted-p':
      return(`
	The <i>p</i>-value of the Fisher Exact test on the
        enrichment of the motif in the foreground vs. the background peptides.
	This value does <b>not</b> accurately represent the statistical 
	significance of the motif, and should be interpreted 
	as a <b>score</b> only.
      `);
    case 'momo-n-tests':
      return(`
	The number of independent tests performed in the search for the motif,
        which is the number of position/residue pairs algorithm motif-x 
	tested for statistical significance.  <b>Note:</b> This field only appears for
	motif-x, and only when the <code>--db-background</code> option is <b>not</b> used.
      `);
    case 'momo-adjusted-p':
      return(`
	The <i>p</i>-value of the Fisher Exact test on the
        enrichment of the motif in the foreground vs. the background 
	peptides, adjusted for the number of independent tests 
	(1 - (1-pvalue)<sup>tests</sup>).  <b>Note 1:</b> This value <b>does</b>
	accurately represent the statistical significance of the motif.
	<b>Note 2:</b> This field only appears for motif-x, and only when the 
	<span style="white-space: nowrap;"><code>--db-background</code></span>
	option is <b>not</b> used.
      `);
    case 'momo-occurrences':
      return(`
	The modified peptides matching the motif.
        <b>Note: </b>For motif-x, the only the the
	occurrences that are not covered by a previous
	motif are shown, so there may be fewer occurrences
	than indicated by <code>FG</code> unless you
	specified the <code>--harvard</code> option.
      `);
    case 'momo-modl-log':
      return(`
	The log of steps taken by the MoDL algorithm in discovering the
	group of motifs.  The log is only shown on the line for the
	first motif in a group output by MoDL.  
	The description length (DL) given the motifs in the group is 
	reported for each step, followed by the motifs in the group.  
	The step (e.g., group of motifs) 
	achieving the lowest description length is preceded by an aserisk ('*').
	<p>
	The steps are followed by three lines stating the
	final step number, the description length of the motifs in 
	the final group, and the decrease in description length
	relative to the empty group (Step 0).
	</p>
      `);

    // AME output fields
    case 'ame-pvalue':
      return(`
	The optimal enrichment <i>p</i>-value of the motif according to the statistical test;
	not adjusted for multiple tests.
      `);
    case 'ame-adj-pvalue':
      return(`
	The optimal enrichment <i>p</i>-value of the motif according to the statistical test,
        adjusted for multiple tests using a Bonferroni correction. ` + extra + `
	If the best <i>p</i>-value is <i>p</i> before adjustment,
        and the number of multiple tests is <i>n</i>, then the adjusted
        <i>p</i>-value is 1 - (1-<i>p</i>)<i><sup>n</sup></i>.
      `);
    case 'ame-evalue':
      return(`
	The expected number of motifs that would be as enriched in the
        (primary) sequences as this one.  The <i>E</i>-value is the adjusted <i>p</i>-value
        multiplied by the number of motifs in the motif file(s).
      `);

    // CentriMo
    case 'centrimo-adj-pvalue':
      return(`
        The statistical significance of the enrichment of the motif, adjusted for multiple tests. ` + extra + `
      `);
    case 'centrimo-evalue':
      var evalue_html = `
        at least one region as enriched for best matches to the motif as the reported region
	(or would have optimal average distance to the sequence center as low as observed, 
	if you used the <code>--cd</code> option).
      `;
      return(`
	The expected number motifs that would have ` + (extra ? extra : evalue_html) + `
	The <i>E</i>-value is the adjusted <i>p</i>-value multiplied by the number of motifs in the
	input files(s).</td> </tr>
      `);
    case 'centrimo-bin-width':
      return(`
        The width (in sequence positions) of the most enriched region (default),
        <b>or</b> two times the average distance between the center of the best site
        and the sequence center if you used the option <code>--cd</code>.
        A best match to the motif is counted as being in the region if
        the center of the motif falls in the region.
      `);
    case 'centrimo-sites-in-bin':
      return(`
	The number of (positive) sequences whose best match to the motif `
	+ (extra ? extra : "falls in the reported region (default) or anywhere in the sequence (if you used the option <code>--cd</code>).") + `
	<br>Note: This number may be less than the number of
	(positive) sequences that have a best match in the region.
	The reason for this is that a sequence may have many matches that score
	equally best.  If <i>n</i> matches have the best score in a sequence, 1/<i>n</i> is added to the
	appropriate bin for each match.</td> </tr>
      `);
    case 'centrimo-mult-tests':
      return(`
	This is the number of multiple tests (<i>n</i>) done for this motif.
	It was used to adjust the <i>p</i>-value of a region for
	multiple tests using the formula:
	<br>&nbsp;&nbsp;
	  <i>p</i>' = 1 - (1-<i>p</i>)<sup><i>n</i></sup> where <i>p</i> is the unadjusted <i>p</i>-value.
	<br>
	The number of multiple tests is the number of regions
	considered times the number of score thresholds considered.
	It depends on the motif length, sequence length, and the type of
	optimizations being done (central enrichment, local enrichment, central distance or
	score optimization).</td> </tr>
      `);

    // GOMo
    case 'gomo-go-term':
      return(`
        The Gene Ontology Consortium term for a specific role or locality.
        Used for annotating genes with their functions.
      `);
    case 'gomo-score':
      return( `
	A score generated as the <a href="https://en.wikipedia.org/wiki/Geometric_mean">
	geometric mean</a> of <a href="https://en.wikipedia.org/wiki/Mann-Whitney_U_test">rank-sum test(s)</a> 
	for the particular Gene Ontology term. The two groups compared by the rank-sum test are scores of genes annotated 
	with the GO term and scores of genes not annotated with the GO term.</td> </tr>
      `);
    case 'gomo-p-value':
      return( `
	An empirically generated <i>p</i>-value for the enrichment of the GO term.<br>
	The null hypothesis is that by shuffling the labels on gene scores, 
	any possible association between the set of genes that a GO term annotates is destroyed. 
	A large number of scores are generated using the null hypothesis and the number of null 
	hypothesis scores that are better than each of the real scores is summed and then divided 
	by the total null hypothesis scores generated to calculate a <i>p</i>-value.</td> </tr>
      `);

    // MAST

    // MCAST
    case 'mcast-cluster-score':
      return( `
        The score that the hidden Markov model created by MCAST assigned to the motif cluster.<br>
        This is the sum of the scores of the individual motif matches in the cluster, plus a
        gap penalty, <i>g</i>, multiplied by the total size of the inter-motif gaps
        in the cluster.  Individual motif match scores are log2(<i>P(s)/p</i>), where <i>s</i> is the 
        <a href="javascript:help_refine('pop_motif_match_score')">log-odds score</a> 
        of the motif match, <i>P(s)</i> is the 
        <a href="javascript:help_refine('pop_motif_match_pvalue')"><i>p</i>-value</a> 
	of the motif match, and <i>p</i> is the user-specified <i>p</i>-value threshold (default: 0.0005).
        <div class="active" id="pop_motif_match_score_act"></div>
        <div class="active" id="pop_motif_match_pvalue_act"></div>
      `);
    case 'mcast-cluster-p-value':
      return( `
	The <i>p</i>-value of the motif cluster score.<br>
	MCAST estimates <i>p</i>-values by fitting an exponential distribution
	to the observed motif cluster scores.
      `);
    case 'mcast-cluster-E-value':
      return( `
        The <i>E</i>-value of the motif cluster score.<br>
	The <i>E</i>-value is an estimate of the <i>number</i> of false positive matches
	with <i>p</i>-values at least as good as this match\'s.  
        MCAST estimates this by multiplying the motif cluster score <i>p</i>-value
        times the (estimated) number of random matches found in the search.
      `);

    // Tomtom
    case 'tomtom-offset':
      return( `
        The offset of the query motif relative to the target motif in the optimal alignment.<br>
 	A positive value indicates the query is shifted right.
      `);
    case 'tomtom-p-value':
      return( `
        The probability that a random motif of the same width as the target would have an 
	optimal alignment with a match score as good or better than the target's.<br>
	Tomtom estimates the <i>p</i>-value using a null model consisting of sampling
	motif columns from all the columns in the set of target motifs.
      `);
    case 'tomtom-E-value':
      return( `
	The expected number of false positives in the matches up to this point.<br>
	Tomtom estimates the <i>E</i>-value by multiplying the <i>p</i>-value by
	the total number of target motifs in all the target databases.
      `);
    case 'tomtom-overlap':
      return( `
	The number of motif columns that overlap in the optimal alignment.
      `);
    case 'tomtom-orientation':
      return( `
	The orientation of the target motif that gave the optimal alignment. ` + extra + `
      `);
  }

} // get_doc_text

//
// Function to replace the innerHTML of element "id" with an HTML paragraph
// containing the text for 'doc_type', which is known to function get_doc_text.
// This function can be used in help pop-ups.
//
function print_doc_para(id, doc_type, extra) {
  html = `<p>` + get_doc_text(doc_type, extra) + `</p>`; 
  document.getElementById(id).insertAdjacentHTML('beforeend', html);
} // print_doc_para
